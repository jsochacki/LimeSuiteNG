#include "TRXLooper.h"

#include "AvgRmsCounter.h"
#include "comms/IDMA.h"
#include "FPGA/FPGA_common.h"
#include "limesuiteng/LMS7002M.h"
#include "limesuiteng/Logger.h"
#include "limesuiteng/StreamMeta.h"
#include "chips/LMS7002M/LMS7002MCSR_Data.h"
#include "protocols/LMSBoards.h"
#include "threadHelper.h"
#include "utilities/DeltaVariable.h"
#include "streaming/DataPacket.h"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <ciso646>
#include <complex>
#include <queue>
#include <fstream>
#include <iostream>
#include <cinttypes>
#include <numeric>

using namespace std::literals::string_literals;

namespace lime {
using namespace LMS7002MCSR_Data;
using namespace std;
using namespace std::chrono;

static constexpr uint16_t defaultSamplesInPkt = 1360;

static constexpr bool showStats{ false };
static constexpr int statsPeriod_ms{ 1000 }; // at 122.88 MHz MIMO, fpga tx pkt counter overflows every 272ms

static constexpr uint32_t k_alignment_tsp_checkpoint_pairs = 8;
static constexpr uint32_t k_alignment_tsp_max_iterations = 128;
static constexpr uint32_t k_alignment_slope_max_iterations = 256;
static constexpr uint32_t k_alignment_quadrature_max_iterations = 128;
static constexpr double k_alignment_quadrature_accept_mean_deg = 45.0;

static_assert(offsetof(FPGA_RxDataPacket, header0) == 0, "unexpected FPGA_RxDataPacket layout");
static_assert(offsetof(FPGA_RxDataPacket, payloadSizeLSB) == 1, "unexpected FPGA_RxDataPacket layout");
static_assert(offsetof(FPGA_RxDataPacket, payloadSizeMSB) == 2, "unexpected FPGA_RxDataPacket layout");
static_assert(offsetof(FPGA_RxDataPacket, reserved) == 3, "unexpected FPGA_RxDataPacket layout");
static_assert(offsetof(FPGA_RxDataPacket, counter) == 8, "unexpected FPGA_RxDataPacket layout");
static_assert(offsetof(FPGA_RxDataPacket, data) == 16, "unexpected FPGA_RxDataPacket layout");
static_assert(sizeof(FPGA_RxDataPacket) == 4096, "unexpected FPGA_RxDataPacket size");

namespace {

std::vector<double> unwrap_phase_degrees(const std::vector<double>& wrapped_phase_degrees)
{
    std::vector<double> unwrapped_phase_degrees = wrapped_phase_degrees;
    if (unwrapped_phase_degrees.empty())
        return unwrapped_phase_degrees;

    for (std::size_t index = 1; index < unwrapped_phase_degrees.size(); ++index)
    {
        double phase_delta_degrees = unwrapped_phase_degrees[index] - unwrapped_phase_degrees[index - 1];
        while (phase_delta_degrees > 180.0)
        {
            unwrapped_phase_degrees[index] -= 360.0;
            phase_delta_degrees -= 360.0;
        }
        while (phase_delta_degrees < -180.0)
        {
            unwrapped_phase_degrees[index] += 360.0;
            phase_delta_degrees += 360.0;
        }
    }
    return unwrapped_phase_degrees;
}

bool linear_fit_phase_vs_bin(const std::vector<int>& bins,
    const std::vector<double>& phase_degrees,
    double* slope_degrees_per_bin,
    double* intercept_degrees,
    double* rms_error_degrees)
{
    if ((bins.size() != phase_degrees.size()) || bins.empty())
        return false;

    double sum_x = 0.0;
    double sum_y = 0.0;
    double sum_xx = 0.0;
    double sum_xy = 0.0;
    const double point_count = static_cast<double>(bins.size());

    for (std::size_t index = 0; index < bins.size(); ++index)
    {
        const double x_value = static_cast<double>(bins[index]);
        const double y_value = phase_degrees[index];
        sum_x += x_value;
        sum_y += y_value;
        sum_xx += x_value * x_value;
        sum_xy += x_value * y_value;
    }

    const double denominator = point_count * sum_xx - sum_x * sum_x;
    if (std::fabs(denominator) < 1.0e-12)
        return false;

    const double fitted_slope = (point_count * sum_xy - sum_x * sum_y) / denominator;
    const double fitted_intercept = (sum_y - fitted_slope * sum_x) / point_count;

    double squared_error_sum = 0.0;
    for (std::size_t index = 0; index < bins.size(); ++index)
    {
        const double x_value = static_cast<double>(bins[index]);
        const double fitted_value = fitted_intercept + fitted_slope * x_value;
        const double error_value = phase_degrees[index] - fitted_value;
        squared_error_sum += error_value * error_value;
    }

    *slope_degrees_per_bin = fitted_slope;
    *intercept_degrees = fitted_intercept;
    *rms_error_degrees = std::sqrt(squared_error_sum / point_count);
    return true;
}

double mean_absolute_value(const std::vector<double>& values)
{
    if (values.empty())
        return 0.0;

    double absolute_sum = 0.0;
    for (double value : values)
        absolute_sum += std::fabs(value);
    return absolute_sum / static_cast<double>(values.size());
}

bool deinterleave_alignment_packet(const StreamConfig& config,
    const FPGA_RxDataPacket& packet,
    std::vector<complex16_t>* channel_a_samples,
    std::vector<complex16_t>* channel_b_samples)
{
    if ((!channel_a_samples) || (!channel_b_samples))
        return false;

    channel_a_samples->assign(2048, complex16_t(0, 0));
    channel_b_samples->assign(2048, complex16_t(0, 0));

    void* destinations[2] = { channel_a_samples->data(), channel_b_samples->data() };
    DataConversion conversion{};
    conversion.srcFormat = config.linkFormat;
    conversion.destFormat = DataFormat::I16;
    conversion.channelCount = 2;

    const uint16_t payload_size_bytes = packet.GetPayloadSize() == 0 ? sizeof(packet.data) : packet.GetPayloadSize();
    const int samples_deinterleaved = Deinterleave(destinations, packet.data, payload_size_bytes, conversion);

    std::fprintf(stderr, "align: samples_deinterleaved=%d payload_bytes=%u\n",
        samples_deinterleaved,
        payload_size_bytes);
    std::fflush(stderr);
    //512 sized packet but header is 16 bits so only 510 samples
    if (samples_deinterleaved < 32)
        return false;

    channel_a_samples->resize(samples_deinterleaved);
    channel_b_samples->resize(samples_deinterleaved);
    return true;
}

static OpStatus WriteAlignmentRegistersToBothChannels(
    LMS7002M* lms,
    uint16_t address_0,
    uint16_t value_0,
    uint16_t address_1,
    uint16_t value_1)
{
    OpStatus status = lms->SetActiveChannel(LMS7002M::Channel::ChA);
    if (status != OpStatus::Success)
        return status;

    status = lms->SPI_write(address_0, value_0, true);
    if (status != OpStatus::Success)
        return status;

    status = lms->SPI_write(address_1, value_1, true);
    if (status != OpStatus::Success)
        return status;

    status = lms->SetActiveChannel(LMS7002M::Channel::ChB);
    if (status != OpStatus::Success)
        return status;

    status = lms->SPI_write(address_0, value_0, true);
    if (status != OpStatus::Success)
        return status;

    status = lms->SPI_write(address_1, value_1, true);
    if (status != OpStatus::Success)
        return status;

    return lms->SetActiveChannel(LMS7002M::Channel::ChA);
}

static double
compute_single_bin_power(const std::vector<complex16_t>& samples,
                         int                             bin,
                         int                             sample_count,
                         int                             dft_length)
{
   const std::complex<double> imaginary_unit(0.0, 1.0);
   const double               pi = std::acos(-1.0);
   std::complex<double>       spectrum(0.0, 0.0);
   for(int sample_index = 0; sample_index < sample_count; ++sample_index)
   {
      const std::complex<double> sample_value(
         static_cast<double>(samples[sample_index].real()),
         static_cast<double>(samples[sample_index].imag()));
      const std::complex<double> phasor
         = std::exp((-2.0 * imaginary_unit * pi * static_cast<double>(bin)
                     * static_cast<double>(sample_index))
                    / static_cast<double>(dft_length));
      spectrum += sample_value * phasor;
   }
   return std::norm(spectrum);
}

static double
max_absolute_value(const std::vector<double>& values)
{
   double maximum_absolute_value = 0.0;
   for(double value : values)
      maximum_absolute_value
         = std::max(maximum_absolute_value, std::fabs(value));
   return maximum_absolute_value;
}

static double
minimum_channel_power(const alignment_bin_result& result)
{
   return std::min(result.power_a, result.power_b);
}

static std::vector<alignment_bin_result>
filter_alignment_bins_by_relative_power(
   const std::vector<alignment_bin_result>& results,
   double                                   keep_within_db)
{
   std::vector<alignment_bin_result> filtered_results;
   double                            strongest_minimum_power = 0.0;
   for(const alignment_bin_result& result : results)
   {
      if(!result.valid) continue;
      strongest_minimum_power
         = std::max(strongest_minimum_power, minimum_channel_power(result));
   }
   if(strongest_minimum_power <= 0.0) return filtered_results;
   const double minimum_allowed_ratio = std::pow(10.0, -keep_within_db / 10.0);
   const double minimum_allowed_power
      = strongest_minimum_power * minimum_allowed_ratio;
   for(const alignment_bin_result& result : results)
   {
      if(!result.valid) continue;
      if(minimum_channel_power(result) >= minimum_allowed_power)
         filtered_results.push_back(result);
   }
   return filtered_results;
}

} // namespace

static struct tm ReadUTC(FPGA* fpga, uint16_t base)
{
    uint32_t addr[3] = { base, base + 1u, base + 2u };
    uint32_t reg[3];
    fpga->ReadRegisters(addr, reg, 3);
    int sec = reg[0] & 0x3F;
    int min = (reg[0] >> 6) & 0x3F;
    int h = reg[1] & 0x1F;
    int d = (reg[1] >> 5) & 0x1F;
    int m = (reg[1] >> 10) & 0xF;
    int y = (reg[2] >> 0) & 0xFFF;

    struct tm tm;
    memset(&tm, 0, sizeof(struct tm));
    tm.tm_sec = sec;
    tm.tm_min = min;
    tm.tm_hour = h;
    tm.tm_mday = d;
    tm.tm_mon = m - 1; // months since January
    tm.tm_year = y - 1900; // years since 1900
    tm.tm_isdst = -1;
    return tm;
}

static int64_t UTC_to_UnixTime(const struct tm& calendarTime)
{
    struct tm datetime;
    memcpy(&datetime, &calendarTime, sizeof(struct tm));
#ifndef __unix__
    time_t unixtime = _mkgmtime(&datetime);
#else
    time_t unixtime = timegm(&datetime);
#endif
    return unixtime;
}

static int ReadySlots(uint32_t writer, uint32_t reader, uint32_t ringSize)
{
    assert(writer < ringSize);
    assert(reader < ringSize);
    if (writer >= reader)
        return writer - reader;
    else
        return ringSize - reader + writer;
}

static constexpr int64_t ts_to_us(int64_t fs, int64_t ts)
{
    int64_t n = (ts / fs);
    int64_t r = (ts % fs);
    return n * 1000000 + ((r * 1000000) / fs);
}

template<class T> static uint32_t indexListToMask(const std::vector<T>& indexes)
{
    uint32_t mask = 0;
    for (T bitIndex : indexes)
        mask |= 1 << bitIndex;
    return mask;
}

static void NegateQChannel(StreamPacket* srcPkt, DataFormat format)
{
    switch (format)
    {
    case DataFormat::I12:
        srcPkt->samples.Scale<lime::complex12_t>(1, -1);
        break;
    case DataFormat::I16:
        srcPkt->samples.Scale<lime::complex16_t>(1, -1);
        break;
    case DataFormat::F32:
        srcPkt->samples.Scale<lime::complex32f_t>(1, -1);
        break;
    default:
        break;
    }
}

/// @brief Constructs a new TRXLooper object.
/// @param rx The DMA communications interface to receive the data from.
/// @param tx The DMA communications interface to send the data to.
/// @param f The FPGA to use in this stream.
/// @param chip The LMS7002M chip to use in this stream.
/// @param moduleIndex The ID of the chip to use.
TRXLooper::TRXLooper(std::shared_ptr<IDMA> rx, std::shared_ptr<IDMA> tx, FPGA* f, LMS7002M* chip, uint8_t moduleIndex)
    : fpga(f)
    , lms(chip)
    , chipId(moduleIndex)
    , mCallback_logMessage(nullptr)
    , mStreamEnabled(false)
    , omitRxPackets(false)
    , startUnixTimeSet(false)
{
    mRxArgs.dma = rx;
    mTxArgs.dma = tx;

    assert(fpga);
    mTimestampOffset = 0;
}

TRXLooper::~TRXLooper()
{
    Stop();
    Teardown();
}

/// @brief Gets the current timestamp of the hardware.
/// @return The current timestamp of the hardware.
uint64_t TRXLooper::GetHardwareTimestamp() const
{
    return mRx.lastTimestamp.load(std::memory_order_relaxed) + mTimestampOffset;
}

/// @brief Sets the hardware timestamp.
/// @param now The current timestamp to set.
/// @return The status of the operation.
OpStatus TRXLooper::SetHardwareTimestamp(const uint64_t now)
{
    mTimestampOffset = now - mRx.lastTimestamp.load(std::memory_order_relaxed);
    return OpStatus::Success;
}

void TRXLooper::Recycle_stream_packets_for_alignment(Stream& stream_state)
{
    StreamPacket* packet_pointer = nullptr;

    if (stream_state.stagingPacket != nullptr)
    {
        stream_state.stagingPacket->Reset();

        if (stream_state.packetsPool)
            stream_state.packetsPool->push(stream_state.stagingPacket, true);
        else
            delete stream_state.stagingPacket;

        stream_state.stagingPacket = nullptr;
    }

    if (stream_state.fifo != nullptr)
    {
        while (stream_state.fifo->pop(&packet_pointer, false))
        {
            if (packet_pointer != nullptr)
            {
                packet_pointer->Reset();

                if (stream_state.packetsPool)
                    stream_state.packetsPool->push(packet_pointer, true);
                else
                    delete packet_pointer;
            }
        }

        stream_state.fifo->clear();
    }
}

OpStatus TRXLooper::Flush_transport_state_for_alignment(void)
{
    OpStatus status = OpStatus::Success;

    if (mStreamEnabled)
    {
        return ReportError(
            OpStatus::Busy, "Flush_transport_state_for_alignment() requires mStreamEnabled == false");
    }

    if (mRx.stage.load(std::memory_order_relaxed) == Stream::ReadyStage::Active)
    {
        return ReportError(OpStatus::Busy, "Flush_transport_state_for_alignment() requires Rx worker idle");
    }

    if (mTx.stage.load(std::memory_order_relaxed) == Stream::ReadyStage::Active)
    {
        return ReportError(OpStatus::Busy, "Flush_transport_state_for_alignment() requires Tx worker idle");
    }

    status = fpga->SelectModule(chipId);
    if (status != OpStatus::Success)
        return status;

    fpga->StopWaveformPlayback();
    fpga->StopStreaming();

    if (mRxArgs.dma)
    {
        status = mRxArgs.dma->Enable(false);
        if (status != OpStatus::Success)
            return status;
    }

    if (mTxArgs.dma)
    {
        status = mTxArgs.dma->Enable(false);
        if (status != OpStatus::Success)
            return status;
    }

    if (mRxArgs.dma)
    {
        for (uint16_t buffer_index = 0; buffer_index < mRxArgs.buffers.size(); ++buffer_index)
            mRxArgs.dma->BufferOwnership(buffer_index, DataTransferDirection::HostToDevice);
    }

    if (mTxArgs.dma)
    {
        for (uint16_t buffer_index = 0; buffer_index < mTxArgs.buffers.size(); ++buffer_index)
            mTxArgs.dma->BufferOwnership(buffer_index, DataTransferDirection::DeviceToHost);
    }

    Recycle_stream_packets_for_alignment(mRx);
    Recycle_stream_packets_for_alignment(mTx);

    mRx.lastTimestamp.store(0, std::memory_order_relaxed);
    mTx.lastTimestamp.store(0, std::memory_order_relaxed);
    mTimestampOffset = 0;

    fpga->ResetPacketCounters(chipId);
    fpga->ResetTimestamp();

    startUnixTime = 0;
    startUnixTimeSet = false;

    mRx.terminate.store(false, std::memory_order_relaxed);
    mTx.terminate.store(false, std::memory_order_relaxed);

    return OpStatus::Success;
}

OpStatus TRXLooper::Discard_initial_rx_dma_transfers_for_alignment(uint32_t number_of_transfers_to_discard, uint8_t irq_period)
{
    OpStatus status = OpStatus::Success;
    IDMA::State dma_state;
    uint64_t last_completed_transfer_count = 0;
    uint64_t submit_request_count = 0;
    uint32_t discarded_transfer_count = 0;
    uint32_t current_buffer_index = 0;
    uint32_t buffer_count = 0;
    uint32_t read_size_bytes = 0;
    bool request_irq = false;

    if (mStreamEnabled)
    {
        return ReportError(
            OpStatus::Busy, "Discard_initial_rx_dma_transfers_for_alignment() requires mStreamEnabled == false");
    }

    if (mRx.stage.load(std::memory_order_relaxed) == Stream::ReadyStage::Active)
    {
        return ReportError(OpStatus::Busy, "Discard_initial_rx_dma_transfers_for_alignment() requires Rx worker idle");
    }

    if (mTx.stage.load(std::memory_order_relaxed) == Stream::ReadyStage::Active)
    {
        return ReportError(OpStatus::Busy, "Discard_initial_rx_dma_transfers_for_alignment() requires Tx worker idle");
    }

    if (mRxArgs.dma == nullptr)
    {
        return ReportError(
            OpStatus::InvalidValue, "Discard_initial_rx_dma_transfers_for_alignment() requires a valid Rx DMA object");
    }

    if (mRxArgs.buffers.empty())
    {
        return ReportError(OpStatus::InvalidValue,
            "Discard_initial_rx_dma_transfers_for_alignment() requires allocated Rx DMA buffers");
    }

    if (mRxArgs.packetSize == 0 || mRxArgs.packetsToBatch == 0)
    {
        return ReportError(
            OpStatus::InvalidValue, "Discard_initial_rx_dma_transfers_for_alignment() requires valid Rx packet sizing");
    }

    if (number_of_transfers_to_discard == 0)
        return OpStatus::Success;

    buffer_count = static_cast<uint32_t>(mRxArgs.buffers.size());
    read_size_bytes = mRxArgs.packetSize * mRxArgs.packetsToBatch;

    dma_state = mRxArgs.dma->GetCounters();
    last_completed_transfer_count = dma_state.transfersCompleted;

    while (discarded_transfer_count < number_of_transfers_to_discard)
    {
        status = mRxArgs.dma->Wait();
        if (status != OpStatus::Success)
            return status;

        dma_state = mRxArgs.dma->GetCounters();

        while (last_completed_transfer_count != dma_state.transfersCompleted &&
            discarded_transfer_count < number_of_transfers_to_discard)
        {
            current_buffer_index = static_cast<uint32_t>(submit_request_count % buffer_count);

            mRxArgs.dma->BufferOwnership(static_cast<uint16_t>(current_buffer_index), DataTransferDirection::DeviceToHost);

            mRxArgs.dma->BufferOwnership(static_cast<uint16_t>(current_buffer_index), DataTransferDirection::HostToDevice);

            request_irq = ((submit_request_count % irq_period) == 0);

            status = mRxArgs.dma->SubmitRequest(
                current_buffer_index, read_size_bytes, DataTransferDirection::DeviceToHost, request_irq);
            if (status != OpStatus::Success)
                return status;

            ++submit_request_count;
            ++discarded_transfer_count;
            ++last_completed_transfer_count;
        }
    }

    return OpStatus::Success;
}

OpStatus TRXLooper::Prepare_rx_transport_for_alignment_capture(void)
{
    OpStatus status = Flush_transport_state_for_alignment();
    if (status != OpStatus::Success)
        return status;

    const uint32_t read_size_bytes = mRxArgs.packetSize;
    constexpr uint8_t irq_period = 1;

    status = mRxArgs.dma->EnableContinuous(true, read_size_bytes, irq_period);
    if (status != OpStatus::Success)
        return status;

    status = fpga->SelectModule(chipId);
    if (status != OpStatus::Success)
    {
        mRxArgs.dma->Enable(false);
        return status;
    }

    fpga->StartStreaming();
    return OpStatus::Success;
}

//OpStatus TRXLooper::Prepare_rx_transport_for_alignment_capture(uint32_t number_of_transfers_to_discard)
//{
//    OpStatus status = Flush_transport_state_for_alignment();
//    uint32_t read_size_bytes = 0;
//    constexpr uint8_t irq_period = 4;
//
//    if (status != OpStatus::Success)
//        return status;
//
//    // More faithful to original limesuite implementation
//    read_size_bytes = sizeof(FPGA_RxDataPacket);
//    //read_size_bytes = mRxArgs.packetSize * mRxArgs.packetsToBatch;
//
//    status = mRxArgs.dma->EnableContinuous(true, read_size_bytes, irq_period);
//    if (status != OpStatus::Success)
//        return status;
//
//    status = fpga->SelectModule(chipId);
//    if (status != OpStatus::Success)
//    {
//        mRxArgs.dma->Enable(false);
//        return status;
//    }
//
//    fpga->StartStreaming();
//
//    status = Discard_initial_rx_dma_transfers_for_alignment(number_of_transfers_to_discard, irq_period);
//    if (status != OpStatus::Success)
//    {
//        fpga->StopStreaming();
//        mRxArgs.dma->Enable(false);
//        return status;
//    }
//
//    return OpStatus::Success;
//}

bool TRXLooper::ShouldAlignRxPhase() const
{
    const auto rx_it = mConfig.channels.find(TRXDir::Rx);
    if (rx_it == mConfig.channels.end())
        return false;

    return mConfig.alignPhase && (rx_it->second.size() == 2);
}

bool TRXLooper::CaptureFreshAlignmentPacket(
    FPGA_RxDataPacket* packet,
    std::chrono::milliseconds timeout)
{
    if (packet == nullptr)
        return false;

    OpStatus status = Flush_transport_state_for_alignment();
    if (status != OpStatus::Success)
        return false;

    status = mRxArgs.dma->Initialize();
    if (status != OpStatus::Success)
        return false;

    const uint16_t buffer_index = 0;

    status = fpga->SelectModule(chipId);
    if (status != OpStatus::Success)
        return false;

    status = mRxArgs.dma->Enable(true);
    if (status != OpStatus::Success)
        return false;

    mRxArgs.dma->BufferOwnership(
        buffer_index,
        DataTransferDirection::HostToDevice);

    const IDMA::State baseline_state = mRxArgs.dma->GetCounters();
    const uint64_t baseline_completed = baseline_state.transfersCompleted;

    status = mRxArgs.dma->SubmitRequest(
        buffer_index,
        mRxArgs.packetSize,
        DataTransferDirection::DeviceToHost,
        true);
    if (status != OpStatus::Success)
    {
        mRxArgs.dma->Enable(false);
        return false;
    }

    fpga->StartStreaming();

    std::fprintf(
        stderr,
        "align: capture baseline completed=%" PRIu64 "\n",
        baseline_completed);
    std::fflush(stderr);

    const auto start_time = std::chrono::steady_clock::now();

    while ((std::chrono::steady_clock::now() - start_time) < timeout)
    {
        const IDMA::State state = mRxArgs.dma->GetCounters();

        if (state.transfersCompleted > baseline_completed)
        {
            std::fprintf(
                stderr,
                "align: capture got completion completed=%" PRIu64 " buffer_index=%u\n",
                state.transfersCompleted,
                static_cast<unsigned>(buffer_index));
            std::fflush(stderr);

            mRxArgs.dma->BufferOwnership(
                buffer_index,
                DataTransferDirection::DeviceToHost);

            std::memset(packet, 0, sizeof(FPGA_RxDataPacket));
            std::memcpy(
                packet,
                mRxArgs.buffers.at(buffer_index),
                mRxArgs.packetSize);

            mRxArgs.dma->BufferOwnership(
                buffer_index,
                DataTransferDirection::HostToDevice);

            fpga->StopStreaming();
            mRxArgs.dma->Enable(false);
            return true;
        }

        std::this_thread::sleep_for(std::chrono::microseconds(50));
    }

    std::fprintf(stderr, "align: capture timeout no completion\n");
    std::fflush(stderr);

    fpga->StopStreaming();
    mRxArgs.dma->Enable(false);
    return false;
}

bool TRXLooper::CaptureFreshAlignmentSamples(
    std::vector<complex16_t>* channel_a_samples,
    std::vector<complex16_t>* channel_b_samples,
    int required_sample_count,
    std::chrono::milliseconds timeout_per_packet)
{
    if ((channel_a_samples == nullptr) || (channel_b_samples == nullptr))
    {
        return false;
    }

    channel_a_samples->clear();
    channel_b_samples->clear();

    OpStatus status = Flush_transport_state_for_alignment();
    if (status != OpStatus::Success)
    {
        return false;
    }

    status = mRxArgs.dma->Initialize();
    if (status != OpStatus::Success)
    {
        return false;
    }

    const uint16_t buffer_index = 0;
    status = fpga->SelectModule(chipId);
    if (status != OpStatus::Success)
    {
        return false;
    }

    status = mRxArgs.dma->Enable(true);
    if (status != OpStatus::Success)
    {
        return false;
    }

    fpga->StartStreaming();

    bool success_flag = false;

    // Continue capturing packets until the vector size reaches the required threshold
    while (static_cast<int>(channel_a_samples->size()) < required_sample_count)
    {
        mRxArgs.dma->BufferOwnership(buffer_index, DataTransferDirection::HostToDevice);

        const IDMA::State baseline_state = mRxArgs.dma->GetCounters();
        const uint64_t baseline_completed = baseline_state.transfersCompleted;

        status = mRxArgs.dma->SubmitRequest(
            buffer_index, 
            mRxArgs.packetSize, 
            DataTransferDirection::DeviceToHost, 
            true);

        if (status != OpStatus::Success)
        {
            break;
        }

        const std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();
        bool packet_received = false;

        while ((std::chrono::steady_clock::now() - start_time) < timeout_per_packet)
        {
            const IDMA::State current_state = mRxArgs.dma->GetCounters();
            if (current_state.transfersCompleted > baseline_completed)
            {
                FPGA_RxDataPacket hardware_packet;
                std::vector<complex16_t> temporary_channel_a;
                std::vector<complex16_t> temporary_channel_b;

                mRxArgs.dma->BufferOwnership(buffer_index, DataTransferDirection::DeviceToHost);
                
                std::memset(&hardware_packet, 0, sizeof(FPGA_RxDataPacket));
                std::memcpy(&hardware_packet, mRxArgs.buffers.at(buffer_index), mRxArgs.packetSize);
                
                mRxArgs.dma->BufferOwnership(buffer_index, DataTransferDirection::HostToDevice);

                // Deinterleave the single packet and append to the aggregate vectors
                if (!deinterleave_alignment_packet(
                        mConfig, 
                        hardware_packet, 
                        &temporary_channel_a, 
                        &temporary_channel_b))
                {
                    packet_received = false;
                    break;
                }

                channel_a_samples->insert(
                    channel_a_samples->end(), 
                    temporary_channel_a.begin(), 
                    temporary_channel_a.end());

                channel_b_samples->insert(
                    channel_b_samples->end(), 
                    temporary_channel_b.begin(), 
                    temporary_channel_b.end());

                packet_received = true;
                break;
            }
            std::this_thread::sleep_for(std::chrono::microseconds(50));
        }

        if (!packet_received)
        {
            break;
        }
    }

    fpga->StopStreaming();
    mRxArgs.dma->Enable(false);

    if (static_cast<int>(channel_a_samples->size()) >= required_sample_count)
    {
        success_flag = true;
    }

    if (success_flag)
    {
        // Truncate to exact required count to maintain DFT alignment if necessary
        channel_a_samples->resize(required_sample_count);
        channel_b_samples->resize(required_sample_count);
    }

    return success_flag;
}

bool TRXLooper::CaptureAlignmentPacket(FPGA_RxDataPacket* packet, std::chrono::milliseconds timeout)
{
    if (packet == nullptr)
        return false;
    const auto start_time = std::chrono::steady_clock::now();
    const auto buffer_count = mRxArgs.buffers.size();
    const IDMA::State initial_state = mRxArgs.dma->GetCounters();
    uint64_t last_completed = initial_state.transfersCompleted;
    std::fprintf(stderr, "align: capture start completed=%" PRIu64 "\n", last_completed);
    std::fflush(stderr);
    while ((std::chrono::steady_clock::now() - start_time) < timeout)
    {
        const IDMA::State state = mRxArgs.dma->GetCounters();
        if (state.transfersCompleted != last_completed)
        {
            const uint64_t completed_index = state.transfersCompleted - 1;
            const uint16_t buffer_index = static_cast<uint16_t>(completed_index % buffer_count);
            std::fprintf(stderr,
                "align: capture got completion completed=%" PRIu64 " buffer_index=%u\n",
                state.transfersCompleted,
                static_cast<unsigned>(buffer_index));
            std::fflush(stderr);
            mRxArgs.dma->BufferOwnership(buffer_index, DataTransferDirection::DeviceToHost);
            std::memset(packet, 0, sizeof(FPGA_RxDataPacket));
            std::memcpy(packet, mRxArgs.buffers.at(buffer_index), mRxArgs.packetSize);
            mRxArgs.dma->BufferOwnership(buffer_index, DataTransferDirection::HostToDevice);
            return true;
        }
        std::this_thread::sleep_for(std::chrono::microseconds(100));
    }
    std::fprintf(stderr, "align: capture timeout no completion\n");
    std::fflush(stderr);
    return false;
}

bool TRXLooper::CheckTSPAligned(const FPGA_RxDataPacket& packet, uint32_t checkpoint_pairs) const
{
    const uint16_t payload_size_bytes = packet.GetPayloadSize() == 0 ? sizeof(packet.data) : packet.GetPayloadSize();
    const uint32_t payload_word_count = payload_size_bytes / sizeof(uint32_t);
    const uint32_t required_word_count = checkpoint_pairs * 2;
    if (payload_word_count < required_word_count)
        return false;

    const uint32_t* payload_words = reinterpret_cast<const uint32_t*>(packet.data);
    for (uint32_t pair_index = 0; pair_index < checkpoint_pairs; ++pair_index)
    {
        const uint32_t left_word = payload_words[2 * pair_index + 0];
        const uint32_t right_word = payload_words[2 * pair_index + 1];
        if (left_word != right_word)
            return false;
    }

    return true;
}

bool TRXLooper::AlignRxTSPRobust(uint32_t checkpoint_pairs)
{
    uint16_t reg0400_a = 0;
    uint16_t reg040c_a = 0;
    uint16_t reg0400_b = 0;
    uint16_t reg040c_b = 0;

    {
        LMS7002M::ChannelScope channel_a_scope(lms, LMS7002M::Channel::ChA);
        reg0400_a = lms->SPI_read(0x0400, true);
        reg040c_a = lms->SPI_read(0x040C, true);
    }
    {
        LMS7002M::ChannelScope channel_b_scope(lms, LMS7002M::Channel::ChB);
        reg0400_b = lms->SPI_read(0x0400, true);
        reg040c_b = lms->SPI_read(0x040C, true);
    }

    {
        const OpStatus write_both_status =
            WriteAlignmentRegistersToBothChannels(lms, 0x0400, 0x8085, 0x040C, 0x01FF);

        if (write_both_status != OpStatus::Success)
        {
            lime::warning("align: failed to write temporary RxTSP alignment registers to both channels");
            return false;
        }
    }


    std::fprintf(stderr, "align: tsp search start\n");
    std::fflush(stderr);

    bool aligned = false;
    for (uint32_t iteration = 0; iteration < k_alignment_tsp_max_iterations; ++iteration)
    {
        {
            LMS7002M::ChannelScope channel_scope(lms, LMS7002M::Channel::ChA);
            lms->SPI_write(0x0020, 0x55FE, true);
            lms->SPI_write(0x0020, 0xFFFD, true);
        }

        //const OpStatus prepare_status = Prepare_rx_transport_for_alignment_capture(2u);
        FPGA_RxDataPacket packet;
        const bool        have_packet
           = CaptureFreshAlignmentPacket(&packet,
                                         std::chrono::milliseconds(150));
        if(!have_packet) continue;

        if (CheckTSPAligned(packet, checkpoint_pairs))
        {
            std::fprintf(stderr, "align: tsp aligned on iteration %u\n", iteration);
            std::fflush(stderr);
            aligned = true;
            break;
        }
    }

    if (!aligned)
    {
        lime::warning("align: tsp search exhausted without success");
    }

    {
        LMS7002M::ChannelScope channel_a_scope(lms, LMS7002M::Channel::ChA);
        lms->SPI_write(0x0400, reg0400_a, true);
        lms->SPI_write(0x040C, reg040c_a, true);
    }
    {
        LMS7002M::ChannelScope channel_b_scope(lms, LMS7002M::Channel::ChB);
        lms->SPI_write(0x0400, reg0400_b, true);
        lms->SPI_write(0x040C, reg040c_b, true);
    }

    return aligned;
}

void TRXLooper::ResetRxIQGeneratorAlignmentState()
{
   const OpStatus channel_status = lms->SetActiveChannel(LMS7002M::Channel::ChA);
   if (channel_status != OpStatus::Success)
   {
       return;
   }

    uint16_t reg20 = lms->SPI_read(0x0020, true);
    uint16_t reg10c = 0;
    uint16_t reg11c = lms->SPI_read(0x011C, true);

    {
        LMS7002M::ChannelScope channel_scope(lms, LMS7002M::Channel::ChA);
        reg10c = lms->SPI_read(0x010C, true);
    }

    lms->SPI_write(0x0020, 0xFFFD, true);
    lms->SPI_write(0x011C, static_cast<uint16_t>(reg11c | 0x0010), true);
    lms->SPI_write(0x0020, 0xFFFF, true);
    lms->SPI_write(0x0124, 0x001F, true);
    lms->SPI_write(0x010C, static_cast<uint16_t>(reg10c | 0x0008), true);
    lms->SPI_write(0x010C, reg10c, true);
    lms->SPI_write(0x0020, 0xFFFD, true);
    lms->SPI_write(0x011C, reg11c, true);
    lms->SPI_write(0x0020, reg20, true);
}

alignment_bin_result TRXLooper::MeasureAlignmentBin(int bin)
{
    alignment_bin_result result;
    result.bin = bin;
    result.valid = false;

    std::vector<complex16_t> channel_a_samples;
    std::vector<complex16_t> channel_b_samples;

    // Aggregating 256 samples (approx 8 packets at 32 samples/packet)
    const int required_sample_count = 256;
    const bool have_samples = CaptureFreshAlignmentSamples(
        &channel_a_samples, 
        &channel_b_samples, 
        required_sample_count, 
        std::chrono::milliseconds(150));

    if (!have_samples)
    {
        return result;
    }

    static constexpr int dft_length = 512;
    const int sample_count = static_cast<int>(channel_a_samples.size());
    const complex64f_t imaginary_unit(0.0, 1.0);
    const double pi_constant = std::acos(-1.0);

    complex64f_t spectrum_a(0.0, 0.0);
    complex64f_t spectrum_b(0.0, 0.0);

    // Phasor: exp(-j * 2 * pi * k * n / N)
    for (int sample_index = 0; sample_index < sample_count; ++sample_index)
    {
        const complex64f_t sample_a(
            static_cast<double>(channel_a_samples[sample_index].real()), 
            static_cast<double>(channel_a_samples[sample_index].imag()));

        const complex64f_t sample_b(
            static_cast<double>(channel_b_samples[sample_index].real()), 
            static_cast<double>(channel_b_samples[sample_index].imag()));

        const double angle = (-2.0 * pi_constant * static_cast<double>(bin) * static_cast<double>(sample_index)) / static_cast<double>(dft_length);
        const complex64f_t phasor = std::exp(imaginary_unit * angle);

        spectrum_a += sample_a * phasor;
        spectrum_b += sample_b * phasor;
    }

    result.power_a = compute_single_bin_power(channel_a_samples, bin, sample_count, dft_length);
    result.power_b = compute_single_bin_power(channel_b_samples, bin, sample_count, dft_length);

    // Calculate phase difference between channels
    result.phase_degrees = (std::arg(spectrum_b) - std::arg(spectrum_a)) * 180.0 / pi_constant;

    // Wrap to [-180, 180]
    while (result.phase_degrees < -180.0) result.phase_degrees += 360.0;
    while (result.phase_degrees > 180.0) result.phase_degrees -= 360.0;

    result.valid = true;
    return result;
}

double
TRXLooper::MeasurePhaseOffsetDeg(int bin, bool* ok)
{
   const alignment_bin_result result = MeasureAlignmentBin(bin);
   if(ok) *ok = result.valid;
   return result.phase_degrees;
}

//double TRXLooper::MeasurePhaseOffsetDeg(int bin, bool* ok)
//{
//    if (ok)
//        *ok = false;
//
//    //const OpStatus prepare_status = Prepare_rx_transport_for_alignment_capture(2u);
//    const OpStatus prepare_status = Prepare_rx_transport_for_alignment_capture(0u);
//    std::fprintf(stderr, "align: prepare_status=%d\n", static_cast<int>(prepare_status));
//    std::fflush(stderr);
//
//    if (prepare_status != OpStatus::Success)
//    {
//        std::fprintf(stderr, "align: prepare failed\n");
//        std::fflush(stderr);
//        return 0.0;
//    }
//
//    FPGA_RxDataPacket packet;
//
//    const bool have_packet = CaptureAlignmentPacket(&packet, std::chrono::milliseconds(50));
//
//   fpga->StopStreaming();
//   mRxArgs.dma->Enable(false);
//
//   std::fprintf(stderr, "align: have_packet=%d\n", have_packet ? 1 : 0);
//   std::fflush(stderr);
//
//   if (!have_packet)
//   {
//       std::fprintf(stderr, "align: no packet captured\n");
//       std::fflush(stderr);
//       return 0.0;
//   }
//
//    std::vector<complex16_t> channel_a_samples;
//    std::vector<complex16_t> channel_b_samples;
//
//    const bool deinterleave_ok =
//    deinterleave_alignment_packet(mConfig, packet, &channel_a_samples, &channel_b_samples);
//
//   std::fprintf(
//       stderr,
//       "align: deinterleave_ok=%d payload_bytes=%u\n",
//       deinterleave_ok ? 1 : 0,
//       packet.GetPayloadSize() == 0
//           ? static_cast<unsigned>(sizeof(packet.data))
//           : packet.GetPayloadSize());
//   std::fflush(stderr);
//
//   if (!deinterleave_ok)
//       return 0.0;
//
//    static constexpr int dft_length = 512;
//    const int sample_count = std::min<int>(512, static_cast<int>(channel_a_samples.size()));
//    const complex64f_t imaginary_unit(0.0, 1.0);
//    const double pi = std::acos(-1.0);
//
//    complex64f_t spectrum_a(0.0, 0.0);
//    complex64f_t spectrum_b(0.0, 0.0);
//    for (int sample_index = 0; sample_index < sample_count; ++sample_index)
//    {
//        const complex64f_t sample_a(channel_a_samples[sample_index].real(), channel_a_samples[sample_index].imag());
//        const complex64f_t sample_b(channel_b_samples[sample_index].real(), channel_b_samples[sample_index].imag());
//        const complex64f_t phasor = std::exp((-2.0 * imaginary_unit * pi * static_cast<double>(bin) *
//                                                 static_cast<double>(sample_index)) /
//            static_cast<double>(dft_length));
//        spectrum_a += sample_a * phasor;
//        spectrum_b += sample_b * phasor;
//    }
//
//   const double power_a
//      = compute_single_bin_power(channel_a_samples, bin, sample_count, dft_length);
//   const double power_b
//      = compute_single_bin_power(channel_b_samples, bin, sample_count, dft_length);
//   double phase_degrees
//      = std::arg(spectrum_b) * 180.0 / pi - std::arg(spectrum_a) * 180.0 / pi;
//   while(phase_degrees < -180.0) phase_degrees += 360.0;
//   while(phase_degrees > 180.0) phase_degrees -= 360.0;
//   std::fprintf(stderr, "align: bin=%d phase_deg=%+.6f power_a=%.3e power_b=%.3e "
//               "payload_bytes=%u samples=%zu\n",
//               bin,
//               phase_degrees,
//               power_a,
//               power_b,
//               packet.GetPayloadSize() == 0
//                  ? static_cast<unsigned>(sizeof(packet.data))
//                  : packet.GetPayloadSize(),
//               channel_a_samples.size());
//   std::fflush(stderr);
//   if(ok) *ok = true;
//   return phase_degrees;
//}

bool TRXLooper::SearchRxPhaseSlopeState(double sample_rate_hz, int decimation_index, const std::vector<int>& bins)
{
    static constexpr double k_alignment_slope_power_keep_within_db = 12.0;
    static constexpr std::size_t k_alignment_slope_min_valid_bins = 4;
    static const double legacy_offsets[] = { 1.15 / 60.0, 1.10 / 40.0, 0.55 / 20.0, 0.20 / 10.0, 0.18 / 5.0 };
    static const double legacy_tolerances[] = { 0.90, 0.45, 0.25, 0.14, 0.06 };

    if (decimation_index < 0 || decimation_index > 4)
        decimation_index = 0;

    const double expected_phase_difference_deg = legacy_offsets[decimation_index] * sample_rate_hz / 1.0e6;
    const double expected_slope_deg_per_bin = -expected_phase_difference_deg / 32.0;
    const double slope_tolerance_deg_per_bin = legacy_tolerances[decimation_index] / 32.0;
    const double residual_rms_tolerance_deg = std::max(2.0, legacy_tolerances[decimation_index] * 8.0);

    std::fprintf(stderr, "align: slope search start\n");
    std::fflush(stderr);

    for (uint32_t iteration = 0; iteration < k_alignment_slope_max_iterations; ++iteration)
    {
        lms->Modify_SPI_Reg_bits(LMS7002MCSR::PD_FDIV_O_CGEN, 1, true);
        lms->Modify_SPI_Reg_bits(LMS7002MCSR::PD_FDIV_O_CGEN, 0, true);

        if (!AlignRxTSPRobust(k_alignment_tsp_checkpoint_pairs))
            continue;


        std::vector<alignment_bin_result> measured_bins;
        measured_bins.reserve(bins.size());
        for(int bin : bins)
        {
           const double tx_frequency_hz
              = 450.0e6 + sample_rate_hz * static_cast<double>(bin) / 512.0;
           lms->SetFrequencySX(TRXDir::Tx, tx_frequency_hz);
           const alignment_bin_result result = MeasureAlignmentBin(bin);
           if(!result.valid)
           {
              measured_bins.clear();
              break;
           }
           measured_bins.push_back(result);
        }
        if(measured_bins.empty()) continue;
        const std::vector<alignment_bin_result> filtered_bins
           = filter_alignment_bins_by_relative_power(
              measured_bins,
              k_alignment_slope_power_keep_within_db);
        if(filtered_bins.size() < k_alignment_slope_min_valid_bins) continue;
        std::vector<int>    filtered_bin_indices;
        std::vector<double> filtered_phase_degrees;
        for(const alignment_bin_result& result : filtered_bins)
        {
           filtered_bin_indices.push_back(result.bin);
           filtered_phase_degrees.push_back(result.phase_degrees);
        }
        const std::vector<double> unwrapped_phase_degrees
           = unwrap_phase_degrees(filtered_phase_degrees);
        double fitted_slope_deg_per_bin = 0.0;
        double fitted_intercept_deg     = 0.0;
        double fitted_rms_error_deg     = 0.0;
        if(!linear_fit_phase_vs_bin(filtered_bin_indices,
                                    unwrapped_phase_degrees,
                                    &fitted_slope_deg_per_bin,
                                    &fitted_intercept_deg,
                                    &fitted_rms_error_deg))
        {
           continue;
        }


        const double slope_error_deg_per_bin =
            std::fabs(fitted_slope_deg_per_bin - expected_slope_deg_per_bin);

        std::fprintf(stderr, 
            "align: slope iter=%u fitted_slope_deg_per_bin=%+.9f expected_slope_deg_per_bin=%+.9f slope_error_deg_per_bin=%.9f rms_error_deg=%.6f\n",
            iteration,
            fitted_slope_deg_per_bin,
            expected_slope_deg_per_bin,
            slope_error_deg_per_bin,
            fitted_rms_error_deg);
           std::fflush(stderr);

        if ((slope_error_deg_per_bin <= slope_tolerance_deg_per_bin) &&
            (fitted_rms_error_deg <= residual_rms_tolerance_deg))
        {
            std::fprintf(stderr, "align: slope search accepted on iteration %u\n", iteration);
           std::fflush(stderr);
            return true;
        }
    }

    lime::warning("align: slope search exhausted without success");
    return false;
}

bool TRXLooper::AlignQuadratureRobust(const std::vector<int>& bins, double accept_abs_mean_phase_deg)
{
   static constexpr double k_alignment_quadrature_power_keep_within_db = 12.0;
   static constexpr std::size_t k_alignment_quadrature_min_valid_bins  = 2;
   const double k_alignment_quadrature_mean_abs_phase_deg = accept_abs_mean_phase_deg;
   static constexpr double k_alignment_quadrature_max_abs_phase_deg    = 20.0;

   auto*                   register_backup = lms->BackupRegisterMap();

   lms->SPI_write(0x0020, 0xFFFF, true);
   lms->SPI_write(0x0113, 0x0046, true);
   lms->SPI_write(0x0118, 0x418C, true);
   lms->SPI_write(0x0100, 0x4039, true);
   lms->SPI_write(0x0101, 0x7801, true);
   lms->SPI_write(0x0108, 0x318C, true);
   lms->SPI_write(0x0082, 0x8001, true);
   lms->SPI_write(0x0200, 0x008D, true);
   lms->SPI_write(0x0208, 0x01FB, true);
   lms->SPI_write(0x0400, 0x8081, true);
   lms->SPI_write(0x040C, 0x01FF, true);
   lms->SPI_write(0x0404, 0x0006, true);

   {
       const OpStatus status = lms->SetActiveChannel(LMS7002M::Channel::ChA);
       if (status != OpStatus::Success)
       {
           if (register_backup)
               lms->RestoreRegisterMap(register_backup);
           return false;
       }
   }

    lms->LoadDC_REG_IQ(TRXDir::Tx, 0x3FFF, 0x3FFF);
    lms->SPI_write(0x0020, 0xFFFE, true);
    lms->SPI_write(0x0105, 0x0006, true);
    lms->SPI_write(0x0100, 0x4038, true);
    lms->SPI_write(0x0113, 0x007F, true);
    lms->SPI_write(0x0119, 0x529B, true);


    uint16_t path_value = lms->Get_SPI_Reg_bits(LMS7002MCSR::SEL_PATH_RFE, true);
    lms->SPI_write(0x010D, path_value == 3 ? 0x018F : path_value == 2 ? 0x0117 : 0x008F, true);
    lms->SPI_write(0x010C, path_value == 2 ? 0x88C5 : 0x88A5, true);
    lms->SPI_write(0x0020, 0xFFFD, true);
    lms->SPI_write(0x0103, path_value == 2 ? 0x0612 : 0x0A12, true);
    path_value = lms->Get_SPI_Reg_bits(LMS7002MCSR::SEL_PATH_RFE, true);
    lms->SPI_write(0x010D, path_value == 3 ? 0x018F : path_value == 2 ? 0x0117 : 0x008F, true);
    lms->SPI_write(0x010C, path_value == 2 ? 0x88C5 : 0x88A5, true);
    lms->SPI_write(0x0119, 0x5293, true);

    const double sample_rate_hz = lms->GetSampleRate(TRXDir::Rx, LMS7002M::Channel::ChA);
    const double rx_frequency_hz = lms->GetFrequencySX(TRXDir::Rx);
    lms->SetFrequencySX(TRXDir::Tx, rx_frequency_hz + sample_rate_hz / 16.0);

    {
        const OpStatus mac_restore_status = lms->SetActiveChannel(LMS7002M::Channel::ChA);
        if (mac_restore_status != OpStatus::Success)
        {
            if (register_backup)
                lms->RestoreRegisterMap(register_backup);
            return false;
        }
    }

    std::fprintf(stderr, "align: forced MAC back to channel A before quadrature search\n");
    std::fprintf(stderr, "align: quadrature search start\n");
    std::fflush(stderr);

    bool aligned = false;
    for (uint32_t iteration = 0; iteration < k_alignment_quadrature_max_iterations; ++iteration)
    {
       std::vector<alignment_bin_result> measured_bins;
       measured_bins.reserve(bins.size());
       for(int bin : bins)
       {
          const alignment_bin_result result = MeasureAlignmentBin(bin);
          if(!result.valid) continue;
          measured_bins.push_back(result);
       }
       const std::vector<alignment_bin_result> filtered_bins
          = filter_alignment_bins_by_relative_power(
             measured_bins,
             k_alignment_quadrature_power_keep_within_db);
       if(filtered_bins.size() < k_alignment_quadrature_min_valid_bins)
       {
          ResetRxIQGeneratorAlignmentState();
          continue;
       }
       std::vector<double> filtered_phase_degrees;
       filtered_phase_degrees.reserve(filtered_bins.size());
       for(const alignment_bin_result& result : filtered_bins)
          filtered_phase_degrees.push_back(result.phase_degrees);
       const std::vector<double> unwrapped_phase_degrees
          = unwrap_phase_degrees(filtered_phase_degrees);
       const double mean_absolute_phase_degrees
          = mean_absolute_value(unwrapped_phase_degrees);
       const double maximum_absolute_phase_degrees
          = max_absolute_value(unwrapped_phase_degrees);
       std::fprintf(stderr,
                    "align: quadrature iter=%u mean_abs_phase_deg=%.6f "
                    "max_abs_phase_deg=%.6f valid_bins=%zu\n",
                    iteration,
                    mean_absolute_phase_degrees,
                    maximum_absolute_phase_degrees,
                    filtered_bins.size());
       std::fflush(stderr);
       if((mean_absolute_phase_degrees
           <= k_alignment_quadrature_mean_abs_phase_deg)
          && (maximum_absolute_phase_degrees
              <= k_alignment_quadrature_max_abs_phase_deg))
       {
          std::fprintf(stderr,
                       "align: quadrature accepted on iteration %u\n",
                       iteration);
          std::fflush(stderr);
          aligned = true;
          break;
       }

    }

    if (!aligned)
    {
        lime::warning("align: quadrature search exhausted without success");
    }

    if (register_backup)
        lms->RestoreRegisterMap(register_backup);

    return aligned;
}

OpStatus TRXLooper::AlignRxPhaseInternal()
{
    if (!ShouldAlignRxPhase())
        return OpStatus::Success;

    const uint16_t mac_backup = lms->SPI_read(0x0020, true);
    auto* register_backup = lms->BackupRegisterMap();

    lms->SPI_write(0x0020, 0xFFFF, true);
    lms->SPI_write(0x010C, 0x88C5, true);
    lms->SPI_write(0x010D, 0x0117, true);
    lms->SPI_write(0x0113, 0x024A, true);
    lms->SPI_write(0x0118, 0x418C, true);
    lms->SPI_write(0x0100, 0x4039, true);
    lms->SPI_write(0x0101, 0x7801, true);
    lms->SPI_write(0x0103, 0x0612, true);
    lms->SPI_write(0x0108, 0x318C, true);
    lms->SPI_write(0x0082, 0x8001, true);
    lms->SPI_write(0x0200, 0x008D, true);
    lms->SPI_write(0x0208, 0x01FB, true);
    lms->SPI_write(0x0400, 0x8081, true);
    lms->SPI_write(0x040C, 0x01FF, true);
    lms->SPI_write(0x0404, 0x0006, true);

    {
        const OpStatus mac_restore_status = lms->SetActiveChannel(LMS7002M::Channel::ChA);
        if (mac_restore_status != OpStatus::Success)
        {
            if (register_backup)
                lms->RestoreRegisterMap(register_backup);

            lms->SPI_write(0x0020, mac_backup, true);
            return mac_restore_status;
        }
    }
    std::fprintf(stderr, "align: forced MAC back to channel A before slope search\n");
    std::fflush(stderr);

    lms->LoadDC_REG_IQ(TRXDir::Tx, 0x3FFF, 0x3FFF);

    const double sample_rate_hz = lms->GetSampleRate(TRXDir::Rx, LMS7002M::Channel::ChA);
    lms->SetFrequencySX(TRXDir::Rx, 450.0e6);

    int decimation_index = lms->Get_SPI_Reg_bits(LMS7002MCSR::HBD_OVR_RXTSP, true);
    if (decimation_index > 4)
        decimation_index = 0;

    std::fprintf(stderr, 
        "align: slope search sample_rate_hz=%.3f decimation_index=%d\n",
        sample_rate_hz,
        decimation_index);
    std::fflush(stderr);

    const std::vector<int> slope_bins = { 16, 24, 32, 40, 48, 56, 64 };
    const bool slope_state_ok = SearchRxPhaseSlopeState(sample_rate_hz, decimation_index, slope_bins);

    if (register_backup)
        lms->RestoreRegisterMap(register_backup);

    lms->SPI_write(0x0020, mac_backup, true);

    if (!slope_state_ok)
    {
        lime::warning("Rx phase alignment failed during slope-state search");
        return OpStatus::Error;
    }

    const std::vector<int> quadrature_bins = { 24, 32, 40 };
    const bool quadrature_ok = AlignQuadratureRobust(quadrature_bins, k_alignment_quadrature_accept_mean_deg);
    if (!quadrature_ok)
    {
        lime::warning("Rx phase alignment failed during quadrature-state search");
        return OpStatus::Error;
    }

    return OpStatus::Success;
}

/// @brief Sets up the stream of this looper.
/// @param cfg The configuration settings to set up the stream with.
/// @return The status of the operation.
OpStatus TRXLooper::Setup(const StreamConfig& cfg)
{
    if (mStreamEnabled)
        return ReportError(OpStatus::Busy, "Samples streaming already running"s);

    // if (cfg.channels.at(lime::TRXDir::Rx).size() > 0 && !mRxArgs.port->IsOpen())
    //     return ReportError(OpStatus::IOFailure, "Rx data port not open"s);
    // if (cfg.channels.at(lime::TRXDir::Tx).size() > 0 && !mTxArgs.port->IsOpen())
    //     return ReportError(OpStatus::IOFailure, "Tx data port not open"s);

    float combinedSampleRate =
        std::max(cfg.channels.at(lime::TRXDir::Tx).size(), cfg.channels.at(lime::TRXDir::Rx).size()) * cfg.hintSampleRate;
    int batchSize = 7; // should be good enough for most cases
    // for high data rates e.g 16bit ADC/DAC 2x2 MIMO @ 122.88Msps = ~1973 MB/s
    // need to batch as many packets as possible into transfer buffers
    if (combinedSampleRate != 0)
    {
        batchSize = combinedSampleRate / 61.44e6;
        batchSize = std::clamp(batchSize, 1, 4);
    }

    if ((cfg.linkFormat != DataFormat::I12) && (cfg.linkFormat != DataFormat::I16))
        return ReportError(OpStatus::InvalidValue, "Unsupported stream link format"s);

    mTx.packetsToBatch = 6;
    mRx.packetsToBatch = 6;

    OpStatus status = fpga->SelectModule(chipId);
    if (status != OpStatus::Success)
        return status;
    fpga->StopStreaming();
    fpga->StopWaveformPlayback();
    fpga->ResetPacketCounters(chipId);
    fpga->ResetTimestamp();

    bool needTx = cfg.channels.at(TRXDir::Tx).size() > 0;
    bool needRx = cfg.channels.at(TRXDir::Rx).size() > 0 || needTx; // always need Rx to know current timestamps, cfg.rxCount > 0;
    omitRxPackets = cfg.channels.at(TRXDir::Rx).size() == 0;

    uint16_t channelEnables = 0;
    channelEnables |= indexListToMask(cfg.channels.at(TRXDir::Rx));
    if (channelEnables & ~0x3)
        return ReportError(OpStatus::InvalidValue, "Invalid Rx channel, only [0,1] channels supported"s);
    channelEnables |= indexListToMask(cfg.channels.at(TRXDir::Tx));
    if (channelEnables & ~0x3)
        return ReportError(OpStatus::InvalidValue, "Invalid Tx channel, only [0,1] channels supported"s);

    mConfig = cfg;

    if (!needTx && !needRx)
        return OpStatus::Success;

    const bool use_trxiqpulse = lms->Get_SPI_Reg_bits(LMS7002MCSR::LML1_TRXIQPULSE);
    const bool sisoddr_on = lms->Get_SPI_Reg_bits(LMS7002MCSR::LML1_SISODDR);
    status = fpga->ConfigureSamplesStream(channelEnables, cfg.linkFormat, sisoddr_on, use_trxiqpulse);
    if (status != OpStatus::Success)
        return status;

    // XTRX has RF switches control bits where the GPS_PPS control should be.
    uint16_t devId = fpga->ReadRegister(0x0000);
    bool hasGPSPPS = devId != LMS_DEV_LIMESDR_XTRX && devId != LMS_DEV_EXTERNAL_SSDR;
    if (hasGPSPPS)
    {
        constexpr uint16_t waitGPS_PPS = 1 << 2;
        int interface_ctrl_000A = fpga->ReadRegister(0x000A);
        interface_ctrl_000A &= ~waitGPS_PPS; // disable by default
        if (cfg.extraConfig.waitPPS)
        {
            interface_ctrl_000A |= waitGPS_PPS;
        }
        fpga->WriteRegister(0x000A, interface_ctrl_000A);
    }

    if (mConfig.timestampType == TimestampType::SAMPLE_TICKS)
    {
        fpga->WriteRegister(0x0280, 0); // samples counting
        ticksPerSample = 1;
    }
    else
    {
        fpga->WriteRegister(0x0280, 1); // PPS, and clock ticks
        uint16_t port1sisoddr = lms->Get_SPI_Reg_bits(LMS7002MCSR::LML1_SISODDR);
        // uint16_t port2sisoddr = lms->Get_SPI_Reg_bits(LMS7002MCSR::LML2_SISODDR);
        ticksPerSample = port1sisoddr ? 1 : 2; // 1 or 2 depending on chip settings
    }

    if (mConfig.extraConfig.waitPPS)
    {
        fpga->WriteRegister(0x0281, 1); // rx start with next PPS
        fpga->WriteRegister(0x0282, 1); // tx delay
    }
    else
    {
        fpga->WriteRegister(0x0281, 0); // rx no delay
        fpga->WriteRegister(0x0282, 0); // tx delay
    }

    RxTeardown();
    if (needRx)
        status = RxSetup();

    if (status != OpStatus::Success)
        return status;

    TxTeardown();
    if (needTx)
        status = TxSetup();

    if (status != OpStatus::Success)
        return status;

    return OpStatus::Success;
}

const StreamConfig& TRXLooper::GetConfig() const
{
    return mConfig;
}

/// @brief Starts the stream of this looper.
OpStatus TRXLooper::Start()
{
    if (mStreamEnabled)
        return OpStatus::Success;

    OpStatus status = fpga->SelectModule(chipId);
    if (status != OpStatus::Success)
        return status;

    // Rx start
    {
        mRx.lastTimestamp.store(0, std::memory_order_relaxed);
        const int32_t readSize = mRxArgs.packetSize * mRxArgs.packetsToBatch;
        constexpr uint8_t irqPeriod{ 4 };
        // Rx DMA has to be enabled before the stream enable, otherwise some data
        // might be lost in the time frame between stream enable and then dma enable.
        mRxArgs.dma->EnableContinuous(true, readSize, irqPeriod);
    }

    if (ShouldAlignRxPhase())
    {
        status = AlignRxPhaseInternal();
        if (status != OpStatus::Success)
        {
            mRxArgs.dma->Enable(false);
            fpga->StopStreaming();
            return status;
        }

        fpga->StopStreaming();
        mRxArgs.dma->Enable(false);
        {
            const int32_t readSize = mRxArgs.packetSize * mRxArgs.packetsToBatch;
            constexpr uint8_t irqPeriod{ 4 };
            mRxArgs.dma->EnableContinuous(true, readSize, irqPeriod);
        }
        fpga->ResetPacketCounters(chipId);
        fpga->ResetTimestamp();
    }
    mRx.terminate.store(false, std::memory_order_relaxed);
    mTx.terminate.store(false, std::memory_order_relaxed);

    fpga->StartStreaming();
    {
        std::lock_guard<std::mutex> lock(streamMutex);
        mStreamEnabled = true;
        streamActive.notify_all();
    }
    return OpStatus::Success;
}

OpStatus TRXLooper::StageStart()
{
    return OpStatus::NotImplemented;
}

/// @brief Stops the stream and cleans up all the memory.
void TRXLooper::Stop()
{
    if (!mStreamEnabled)
        return;
    lime::debug("TRXLooper::Stop()");
    mStreamEnabled = false;

    // wait for loop ends
    if (mRx.stage.load(std::memory_order_relaxed) != Stream::ReadyStage::Disabled)
    {
        mRx.terminate.store(true, std::memory_order_relaxed);
        lime::debug("TRXLooper: wait for Rx loop end.");
        {
            std::unique_lock lck{ mRx.mutex };
            while (mRx.stage.load(std::memory_order_relaxed) == Stream::ReadyStage::Active)
                mRx.cv.wait(lck);
        }
        mRxArgs.dma->Enable(false);

        if (mCallback_logMessage)
        {
            char msg[256];
            std::snprintf(msg, sizeof(msg), "Rx%i stop: packetsIn: %" PRIi64, chipId, mRx.stats.packets);
            mCallback_logMessage(LogLevel::Verbose, msg);
        }
    }

    // wait for loop ends
    if (mTx.stage.load(std::memory_order_relaxed) != Stream::ReadyStage::Disabled)
    {
        mTx.terminate.store(true, std::memory_order_relaxed);
        lime::debug("TRXLooper: wait for Tx loop end."s);
        {
            std::unique_lock lck{ mTx.mutex };
            while (mTx.stage.load(std::memory_order_relaxed) == Stream::ReadyStage::Active)
                mTx.cv.wait(lck);
        }
        mTxArgs.dma->Enable(false);

        uint32_t fpgaTxPktIngressCount;
        uint32_t fpgaTxPktDropCounter;
        fpga->ReadTxPacketCounters(chipId, &fpgaTxPktIngressCount, &fpgaTxPktDropCounter);
        if (mCallback_logMessage)
        {
            char msg[512];
            std::snprintf(msg,
                sizeof(msg),
                "Tx%i stop: host sent packets: %" PRIi64 " (0x%08" PRIX64 "), FPGA packet ingresed: %i (0x%08X), diff: %" PRIi64
                ", Tx packet dropped: %i",
                chipId,
                mTx.stats.packets,
                mTx.stats.packets,
                fpgaTxPktIngressCount,
                fpgaTxPktIngressCount,
                (mTx.stats.packets & 0xFFFFFFFF) - fpgaTxPktIngressCount,
                fpgaTxPktDropCounter);
            mCallback_logMessage(LogLevel::Verbose, msg);
        }
    }

    // Disable FPGA streaming only after data transfer threads finish work.
    // Becase stream disable halts DMA, and threads could get stuck waiting for interrupt
    // of the next data batch.
    fpga->StopStreaming();

    if (mRx.stagingPacket != nullptr)
    {
        mRx.packetsPool->push(mRx.stagingPacket, true);
        mRx.stagingPacket = nullptr;
    }
    if (mRx.fifo)
    {
        while (mRx.fifo->pop(&mRx.stagingPacket, false))
            mRx.packetsPool->push(mRx.stagingPacket, true);
        mRx.fifo->clear();
        mRx.stagingPacket = nullptr;
    }
    if (mTx.stagingPacket != nullptr)
    {
        mTx.packetsPool->push(mTx.stagingPacket, true);
        mTx.stagingPacket = nullptr;
    }
    if (mTx.fifo)
    {
        while (mTx.fifo->pop(&mTx.stagingPacket, false))
            mTx.packetsPool->push(mTx.stagingPacket, true);
        mTx.fifo->clear();
        mTx.stagingPacket = nullptr;
    }

    mRx.lastTimestamp.store(0, std::memory_order_relaxed);
    fpga->ResetPacketCounters(chipId);
    fpga->ResetTimestamp();
    startUnixTimeSet = false;
    startUnixTime = 0;
}

/// @brief Stops all the running streams and clears up the memory.
void TRXLooper::Teardown()
{
    RxTeardown();
    TxTeardown();
}

OpStatus TRXLooper::RxSetup()
{
    OpStatus status = mRxArgs.dma->Initialize();
    if (status != OpStatus::Success)
        return status;

    mRx.terminate.store(false, std::memory_order_relaxed);

    mRx.lastTimestamp.store(0, std::memory_order_relaxed);
    const bool usePoll = mConfig.extraConfig.usePoll;
    const int chCount = std::max(mConfig.channels.at(lime::TRXDir::Rx).size(), mConfig.channels.at(lime::TRXDir::Tx).size());
    assert(chCount > 0);
    const int sampleSize = (mConfig.linkFormat == DataFormat::I16 ? 4 : 3); // sizeof IQ pair

    constexpr std::size_t headerSize{ sizeof(StreamHeader) };

    const GatewareFeatures gw = fpga->GetFeatures();
    uint32_t packetSize = 4096;
    if (gw.hasConfigurableStreamPacketSize)
    {
        int requestSamplesInPkt = 256 / chCount;
        if (mConfig.extraConfig.rx.samplesInPacket > 0)
        {
            requestSamplesInPkt = mConfig.extraConfig.rx.samplesInPacket;
        }
        packetSize = fpga->SetUpVariableRxSize(requestSamplesInPkt, sampleSize, chCount, chipId);
        mRx.samplesInPkt = (packetSize - headerSize) / (sampleSize * chCount);
    }
    else
    {
        mRx.samplesInPkt = (packetSize - headerSize) / (sampleSize * chCount);
    }

    const auto dmaChunks{ mRxArgs.dma->GetBuffers() };
    const auto dmaBufferSize = dmaChunks.front().size;

    if (mConfig.hintSampleRate == 0)
    {
        uint8_t samplerateChannel = 0;
        // if no Rx channels are configured for streaming use channel 0 as reference sample rate
        if (!mConfig.channels.at(TRXDir::Rx).empty())
            samplerateChannel = mConfig.channels.at(TRXDir::Rx).at(0);

        mConfig.hintSampleRate =
            lms->GetSampleRate(TRXDir::Rx, samplerateChannel == 0 ? LMS7002M::Channel::ChA : LMS7002M::Channel::ChB);
    }

    // aim batch size to desired data output period, ~100us should be good enough
    if (mConfig.hintSampleRate > 0)
        mRx.packetsToBatch = std::floor((0.0001 * mConfig.hintSampleRate) / mRx.samplesInPkt);

    if (mConfig.extraConfig.rx.packetsInBatch != 0)
        mRx.packetsToBatch = mConfig.extraConfig.rx.packetsInBatch;

    mRx.packetsToBatch = std::clamp<uint32_t>(mRx.packetsToBatch, 1u, dmaBufferSize / packetSize);

    float bufferTimeDuration;
    if (mConfig.hintSampleRate)
        bufferTimeDuration = mRx.samplesInPkt * mRx.packetsToBatch / mConfig.hintSampleRate;
    else
        bufferTimeDuration = 0;

    std::vector<uint8_t*> dmaBuffers(dmaChunks.size());
    for (uint32_t i = 0; i < dmaChunks.size(); ++i)
    {
        dmaBuffers[i] = dmaChunks[i].buffer;
    }

    mRxArgs.buffers = std::move(dmaBuffers);
    mRxArgs.bufferSize = dmaBufferSize;
    mRxArgs.packetSize = packetSize;
    mRxArgs.packetsToBatch = mRx.packetsToBatch;
    mRxArgs.samplesInPacket = mRx.samplesInPkt;

    assert(mRxArgs.bufferSize > 0);
    assert(mRxArgs.packetSize > 0);
    assert(mRxArgs.packetsToBatch > 0);
    assert(mRxArgs.samplesInPacket > 0);

    const int packetsInFIFO = 0.25 * mConfig.hintSampleRate / mRx.samplesInPkt; // buffer 0.25 second of data
    mRx.packetsPool = std::make_unique<PacketsFIFO<StreamPacket*>>(packetsInFIFO);
    const uint32_t userSampleSize = mConfig.format == DataFormat::F32 ? sizeof(lime::complex32f_t) : sizeof(lime::complex16_t);
    for (uint32_t i = 0; i < mRx.packetsPool->max_size(); ++i)
        mRx.packetsPool->push(new StreamPacket(mRx.samplesInPkt, chCount, userSampleSize));
    mRx.fifo = std::make_unique<PacketsFIFO<StreamPacket*>>(packetsInFIFO);

    char msg[256];
    std::snprintf(msg,
        sizeof(msg),
        "%s Rx%i Setup: usePoll:%i rxSamplesInPkt:%i rxPacketsInBatch:%i, DMA_ReadSize:%i, link:%s, batchSizeInTime:%gus FS:%f, "
        "FIFO=%i*%i\n",
        mRxArgs.dma->GetName().c_str(),
        chipId,
        usePoll ? 1 : 0,
        mRx.samplesInPkt,
        mRx.packetsToBatch,
        mRx.packetsToBatch * packetSize,
        (mConfig.linkFormat == DataFormat::I12 ? "I12" : "I16"),
        bufferTimeDuration * 1e6,
        mConfig.hintSampleRate,
        packetsInFIFO,
        mRx.samplesInPkt);
    if (showStats)
        printf("%s", msg);
    if (mCallback_logMessage)
        mCallback_logMessage(LogLevel::Verbose, msg);

    // Don't just use REALTIME scheduling, or at least be cautious with it.
    // if the thread blocks for too long, Linux can trigger RT throttling
    // which can cause unexpected data packet losses and timing issues.
    // Also need to set policy to default here, because if host process is running
    // with REALTIME policy, these threads would inherit it and exhibit mentioned
    // issues.
    const auto schedulingPolicy = ThreadPolicy::REALTIME;
    mRx.terminate.store(false, std::memory_order_relaxed);
    mRx.terminateWorker.store(false, std::memory_order_relaxed);

    auto RxLoopFunction = std::bind(&TRXLooper::RxWorkLoop, this);
    mRx.thread = std::thread(RxLoopFunction);
    SetOSThreadPriority(ThreadPriority::HIGHEST, schedulingPolicy, &mRx.thread);
#ifdef __linux__
    char threadName[16]; // limited to 16 chars, including null byte.
    snprintf(threadName, sizeof(threadName), "lime:Rx%i", chipId);
    pthread_setname_np(mRx.thread.native_handle(), threadName);
#endif

    // wait for Rx thread to be ready
    lime::debug("RxSetup wait for Rx worker thread."s);
    {
        std::unique_lock lck{ mRx.mutex };
        while (mRx.stage.load(std::memory_order_relaxed) < Stream::ReadyStage::WorkerReady)
            mRx.cv.wait(lck);
    }

    return status;
}

struct DMATransactionCounter {
    uint64_t requests{ 0 };
    uint64_t completed{ 0 };
};

void TRXLooper::RxWorkLoop()
{
    lime::debug("Rx worker thread ready.");
    // signal that thread is ready for work
    {
        std::unique_lock lck{ mRx.mutex };

        // signal that thread is ready for work
        mRx.stage.store(Stream::ReadyStage::WorkerReady, std::memory_order_relaxed);
        mRx.cv.notify_all();
    }

    while (!mRx.terminateWorker.load(std::memory_order_relaxed))
    {
        // thread ready for work, just wait for stream enable
        {
            std::unique_lock lk{ streamMutex };
            while (!mStreamEnabled && !mRx.terminateWorker.load(std::memory_order_relaxed))
                streamActive.wait_for(lk, std::chrono::milliseconds(100));
        }
        if (!mStreamEnabled)
            continue;

        mRx.stage.store(Stream::ReadyStage::Active, std::memory_order_relaxed);
        ReceivePacketsLoop();

        std::unique_lock lck{ mRx.mutex };
        mRx.stage.store(Stream::ReadyStage::WorkerReady, std::memory_order_relaxed);
        mRx.cv.notify_all();
    }
    mRx.stage.store(Stream::ReadyStage::Disabled, std::memory_order_relaxed);
    lime::debug("Rx worker thread shutdown.");
}

static Timespec ExtractPacketTimestamp(const StreamConfig& config, const FPGA_RxDataPacket* fpgapacket, int clockTicksPerSample)
{
    switch (config.timestampType)
    {
    case TimestampType::SAMPLE_TICKS:
        return Timespec(0, fpgapacket->counter, config.hintSampleRate);
        break;
    case TimestampType::REALTIME_SECONDS:
    case TimestampType::UNIX_EPOCH: {
        uint32_t clockCount = fpgapacket->counter & 0xFFFFFFFF;
        uint32_t PPScount = (fpgapacket->counter >> 32) & 0xFFFFFFFF;
        uint64_t ticks = clockCount; // depending on interface configuration there might be 2 or 1 tick per sample
        uint64_t seconds = PPScount;
        return Timespec(seconds, ticks, clockTicksPerSample * config.hintSampleRate);
    }
    break;
    }
    return Timespec();
}

static int32_t ExtractPacketSamples(
    const StreamConfig& config, const TRXLooper::TransferArgs& args, StreamPacket* userPkt, const FPGA_RxDataPacket* fpgapacket)
{
    DataConversion conversion{};
    conversion.srcFormat = config.linkFormat;
    conversion.destFormat = config.format;
    conversion.channelCount = std::max(config.channels.at(lime::TRXDir::Tx).size(), config.channels.at(lime::TRXDir::Rx).size());

    assert(userPkt);
    assert(userPkt->samples.isFull() == false);

    const size_t payloadSize{ args.packetSize - sizeof(StreamHeader) };
    const int samplesProduced = Deinterleave(userPkt->samples.back(), fpgapacket->data, payloadSize, conversion);
    userPkt->samples.SetSize(userPkt->samples.size() + samplesProduced);
    return samplesProduced;
}

static std::string timespec_to_utc_string(struct timespec* ts)
{
    const int nanosecond_offset = 21;
    char buf[64];
    struct tm tm;
#ifdef __unix__
    gmtime_r(&ts->tv_sec, &tm);
#else
    gmtime_s(&tm, &ts->tv_sec);
#endif
    strftime(buf, nanosecond_offset, "%Y-%m-%dT%H:%M:%S.", &tm);
    sprintf(buf + nanosecond_offset - 1, "%09luZ", ts->tv_nsec);
    return std::string(buf);
}

static std::string TimestampToString(Timespec timestamp, TimestampType type)
{
    char timestampstr[256];
    switch (type)
    {
    case TimestampType::SAMPLE_TICKS:
        snprintf(timestampstr, sizeof(timestampstr), "%" PRIu64, timestamp.GetTicks());
        break;
    case TimestampType::REALTIME_SECONDS:
        snprintf(timestampstr, sizeof(timestampstr), "%.9fs", timestamp.GetRealSeconds());
        break;
    case TimestampType::UNIX_EPOCH: {
        struct timespec ts;
        ts.tv_sec = timestamp.GetSeconds();
        ts.tv_nsec = timestamp.GetFracSeconds() * 1e9;
        std::string utctimestamp = timespec_to_utc_string(&ts);
        snprintf(timestampstr, sizeof(timestampstr), "%s", utctimestamp.c_str());
        break;
    }
    }
    return std::string(timestampstr);
}

/** @brief Function dedicated for receiving data samples from board */
void TRXLooper::ReceivePacketsLoop()
{
    lime::debug("Rx receive loop start.");

    const int32_t bufferCount = mRxArgs.buffers.size();
    const int32_t readSize = mRxArgs.packetSize * mRxArgs.packetsToBatch;
    const int32_t packetSize = mRxArgs.packetSize;
    const std::vector<uint8_t*>& dmaBuffers{ mRxArgs.buffers };
    StreamStats& stats = mRx.stats;
    auto& fifo = mRx.fifo;

    DeltaVariable<int32_t> overrun(0);
    DeltaVariable<int32_t> loss(0);

    constexpr uint8_t irqPeriod{ 4 };

    auto t1{ std::chrono::steady_clock::now() };
    auto t2 = t1;

    int32_t Bps = 0;
    StreamPacket* userPkt = nullptr;

    uint32_t lastHwIndex{ 0 };
    DMATransactionCounter counters;

    assert(mRx.stagingPacket == nullptr); // should be clean start
    assert(fifo->empty());

    bool getStartTime = mConfig.timestampType != TimestampType::SAMPLE_TICKS;
    startUnixTime = 0;

    Timespec lastPacketTS;
    Timespec expectedTimestamp;
    Timespec fpgaFrontEndDelay;
    expectedTimestamp.SetTickRate(ticksPerSample * mConfig.hintSampleRate);

    while (mRx.terminate.load(std::memory_order_relaxed) == false)
    {
        IDMA::State dma{ mRxArgs.dma->GetCounters() };
        int64_t counterDiff = ReadySlots(dma.transfersCompleted, lastHwIndex, 65536);
        lastHwIndex = dma.transfersCompleted;
        counters.completed += counterDiff;

        if (counterDiff > 0)
        {
            const int bytesTransferred = counterDiff * readSize;
            Bps += bytesTransferred;
            stats.bytesTransferred += bytesTransferred;
        }

        // print stats
        t2 = std::chrono::steady_clock::now();
        const auto timePeriod{ std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() };
        if (timePeriod >= statsPeriod_ms)
        {
            t1 = t2;
            double dataRateBps = 1000.0 * Bps / timePeriod;
            stats.dataRate_Bps = dataRateBps;
            char msg[512];
            std::snprintf(msg,
                sizeof(msg) - 1,
                "%s Rx%i: %3.3f MB/s | TS:%s pkt:%" PRIi64 " o:%i(%+i) l:%i(%+i) dma:%" PRIu64 "/%" PRIu64 "(+%" PRIu64
                ") swFIFO:%" PRIuPTR,
                mRxArgs.dma->GetName().c_str(),
                chipId,
                stats.dataRate_Bps / 1e6,
                TimestampToString(lastPacketTS, mConfig.timestampType).c_str(),
                stats.packets,
                overrun.value(),
                overrun.delta(),
                loss.value(),
                loss.delta(),
                counters.requests,
                counters.completed,
                counters.completed - counters.requests,
                fifo->size());
            if (showStats)
                printf("%s\n", msg);
            if (mCallback_logMessage)
            {
                bool showAsWarning = overrun.delta() || loss.delta();
                LogLevel level = showAsWarning ? LogLevel::Warning : LogLevel::Debug;
                mCallback_logMessage(level, msg);
            }
            overrun.checkpoint();
            loss.checkpoint();
            Bps = 0;
        }

        if (counters.completed - counters.requests == 0)
        {
            if (mConfig.extraConfig.usePoll)
                mRxArgs.dma->Wait();
            else
                std::this_thread::yield();
            continue;
        }

        const uint64_t currentBufferIndex{ counters.requests % bufferCount };
        mRxArgs.dma->BufferOwnership(currentBufferIndex, DataTransferDirection::DeviceToHost);
        const uint8_t* buffer{ dmaBuffers.at(currentBufferIndex) };

        if (getStartTime)
        {
            std::lock_guard<std::mutex> lock(startTimeMutex);
            startUnixTime = UTC_to_UnixTime(ReadUTC(fpga, 0x0283));
            ++startUnixTime; // last GNSS message is from previous PPS, so add 1 second
            startUnixTimeSet = true;
            // stream starts with the PPS signal, but there is delay until the samples are put into packet by the hardware
            // rewind timstamps by the amount of first packet clockCounter
            if (mConfig.extraConfig.waitPPS)
            {
                const FPGA_RxDataPacket* hardwarePkt = reinterpret_cast<const FPGA_RxDataPacket*>(buffer);
                fpgaFrontEndDelay = ExtractPacketTimestamp(mConfig, hardwarePkt, ticksPerSample);
                expectedTimestamp = expectedTimestamp + fpgaFrontEndDelay;
            }
            getStartTime = false;
            t1 = std::chrono::steady_clock::now();
        }
        startTimeIsSet.notify_all();

        bool reportProblems = false;
        const int srcPktCount = mRxArgs.packetsToBatch;
        for (int i = 0; i < srcPktCount; ++i)
        {
            const FPGA_RxDataPacket* hardwarePkt = reinterpret_cast<const FPGA_RxDataPacket*>(&buffer[packetSize * i]);
            Timespec hwts = ExtractPacketTimestamp(mConfig, hardwarePkt, ticksPerSample);

            // clock counter can drift, and gets reset with PPS, creating small discontinuity in expected and received counter
            Timespec diff = (hwts - expectedTimestamp);
            const Timespec samplePeriod(0, 1.0 / mConfig.hintSampleRate);
            if (abs(diff) > samplePeriod)
            {
                int64_t fpgaTicks = hardwarePkt->counter & 0xFFFFFFFF;
                if (mConfig.timestampType != TimestampType::SAMPLE_TICKS)
                {
                    lime::debug("Expected PPS=%li, CLK=%li, got FPGA PPS=%li, CLK=%li, diff=%lins",
                        expectedTimestamp.GetSeconds(),
                        static_cast<int64_t>(expectedTimestamp.GetFracSeconds() * ticksPerSample * mConfig.hintSampleRate),
                        hardwarePkt->counter >> 32,
                        fpgaTicks,
                        static_cast<int64_t>((diff.GetSeconds() + diff.GetFracSeconds()) * 1e9));
                }
                lime::debug("Loss: pkt:%li exp: %016lx, got: %016lx, diff: %+li, timeDiff:%+lins",
                    stats.packets + i,
                    expectedTimestamp.GetTicks(),
                    hwts.GetTicks(),
                    hwts.GetTicks() - expectedTimestamp.GetTicks(),
                    static_cast<int64_t>((diff.GetSeconds() + diff.GetFracSeconds()) * 1e9));
                ++stats.loss;
                loss.add(1);
                reportProblems = true;
            }
            if (hardwarePkt->txWasDropped())
            {
                ++mTx.stats.loss;
            }
            expectedTimestamp = hwts;
            expectedTimestamp.AddTicks(mRxArgs.samplesInPacket * ticksPerSample);

            lastPacketTS = hwts;
            if (mConfig.timestampType == TimestampType::UNIX_EPOCH)
                lastPacketTS = lastPacketTS + Timespec(startUnixTime);

            if (omitRxPackets)
                continue;

            if (userPkt == nullptr)
            {
                if (!mRx.packetsPool->pop(&userPkt, false) || userPkt == nullptr)
                {
                    ++stats.overrun;
                    overrun.add(1);
                    reportProblems = true;
                    continue;
                }
                userPkt->Reset();
                userPkt->meta.timestamp = hwts;
                userPkt->meta.useTimestamp = true;
                userPkt->meta.flush = false;
            }

            ExtractPacketSamples(mConfig, mRxArgs, userPkt, hardwarePkt);

            if (mConfig.extraConfig.negateQ)
                NegateQChannel(userPkt, mConfig.format);

            if (mConfig.timestampType == TimestampType::UNIX_EPOCH)
            {
                userPkt->meta.timestamp = userPkt->meta.timestamp + Timespec(startUnixTime);
                userPkt->meta.timestamp = userPkt->meta.timestamp - fpgaFrontEndDelay;
            }

            if (fifo->push(userPkt, false))
                userPkt = nullptr;
            else
            {
                ++stats.overrun;
                overrun.add(1);
                userPkt->Reset();
                reportProblems = true;
            }
        }

        stats.packets += srcPktCount;
        stats.timestamp = expectedTimestamp.GetTicks();
        mRx.lastTimestamp.store(expectedTimestamp.GetTicks(), std::memory_order_relaxed);

        mRxArgs.dma->BufferOwnership(currentBufferIndex, DataTransferDirection::HostToDevice);
        bool requestIRQ = (counters.requests % irqPeriod) == 0;
        ++counters.requests;
        mRxArgs.dma->SubmitRequest(currentBufferIndex, readSize, DataTransferDirection::DeviceToHost, requestIRQ);

        // one callback for the entire batch
        if (reportProblems && mConfig.statusCallback)
            mConfig.statusCallback(false, &stats, mConfig.userData);
        std::this_thread::yield();
    }
    lime::debug("Rx receive loop end.");
}

void TRXLooper::RxTeardown()
{
    if (mRx.stage.load(std::memory_order_relaxed) != Stream::ReadyStage::Disabled)
    {
        lime::debug("RxTeardown wait for Rx worker shutdown.");
        mRx.terminateWorker.store(true, std::memory_order_relaxed);
        {
            std::unique_lock lck{ streamMutex };
            mRx.terminateWorker.store(true, std::memory_order_relaxed);
            streamActive.notify_all();
        }

        mRx.terminate.store(true, std::memory_order_relaxed);
        try
        {
            if (mRx.thread.joinable())
                mRx.thread.join();
        } catch (...)
        {
            lime::error("Failed to join TRXLooper Rx thread"s);
        }
    }

    if (mRx.stagingPacket)
    {
        delete mRx.stagingPacket;
        mRx.stagingPacket = nullptr;
    }

    if (mRx.packetsPool)
    {
        while (mRx.packetsPool->pop(&mRx.stagingPacket, false))
            delete mRx.stagingPacket;
        mRx.stagingPacket = nullptr;
    }

    delete mRx.fifo.release();
    delete mRx.packetsPool.release();
}

template<class T>
uint32_t TRXLooper::StreamRxTemplate(T* const* dest, uint32_t count, StreamRxMeta* meta, chrono::microseconds timeout)
{
    bool timestampSet = false;
    uint32_t samplesProduced = 0;
    const bool useChannelB = mConfig.channels.at(TRXDir::Rx).size() > 1;

    bool firstIteration = true;

    assert(dest);
    assert(dest[0]);
    if (useChannelB)
        assert(dest[1]);

    auto start = chrono::high_resolution_clock::now();
    while (samplesProduced < count)
    {
        if (!mRx.stagingPacket && !mRx.fifo->pop(&mRx.stagingPacket, firstIteration, timeout))
        {
            lime::error("No samples or timeout"s);
            return samplesProduced;
        }

        if (!timestampSet && meta)
        {
            meta->timestamp = mRx.stagingPacket->meta.timestamp;
            meta->hasTimestamp = true;
            timestampSet = true;
        }

        uint32_t expectedCount = count - samplesProduced;
        const uint32_t samplesToCopy = std::min(expectedCount, mRx.stagingPacket->samples.size());

        T* const* src = reinterpret_cast<T* const*>(mRx.stagingPacket->samples.front());

        std::memcpy(&dest[0][samplesProduced], src[0], samplesToCopy * sizeof(T));

        if (useChannelB)
        {
            assert(dest[1]);
            std::memcpy(&dest[1][samplesProduced], src[1], samplesToCopy * sizeof(T));
        }

        mRx.stagingPacket->samples.pop(samplesToCopy);
        mRx.stagingPacket->meta.timestamp.AddTicks(samplesToCopy * ticksPerSample);

        samplesProduced += samplesToCopy;

        if (mRx.stagingPacket->samples.empty())
        {
            mRx.packetsPool->push(mRx.stagingPacket);
            mRx.stagingPacket = nullptr;
        }

        auto duration = chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - start);
        if (duration > timeout)
            return samplesProduced;
    }

    return samplesProduced;
}

/// @brief Receives samples from this specific stream.
/// @param samples The buffer to put the received samples in.
/// @param count The amount of samples to receive.
/// @param meta The metadata of the packets of the stream.
/// @return The amount of samples received.
uint32_t TRXLooper::StreamRx(
    lime::complex32f_t* const* samples, uint32_t count, StreamMeta* meta, std::chrono::microseconds timeout)
{
    StreamRxMeta rxmeta;
    uint32_t samplesRead = StreamRxTemplate(samples, count, &rxmeta, timeout);
    if (meta)
        meta->timestamp = rxmeta.timestamp.GetTicks();
    return samplesRead;
}

/// @copydoc TRXLooper::StreamRx()
uint32_t TRXLooper::StreamRx(lime::complex16_t* const* samples, uint32_t count, StreamMeta* meta, std::chrono::microseconds timeout)
{
    StreamRxMeta rxmeta;
    uint32_t samplesRead = StreamRxTemplate(samples, count, &rxmeta, timeout);
    if (meta)
        meta->timestamp = rxmeta.timestamp.GetTicks();
    return samplesRead;
}

/// @copydoc TRXLooper::StreamRx()
uint32_t TRXLooper::StreamRx(lime::complex12_t* const* samples, uint32_t count, StreamMeta* meta, std::chrono::microseconds timeout)
{
    StreamRxMeta rxmeta;
    uint32_t samplesRead = StreamRxTemplate(samples, count, &rxmeta, timeout);
    if (meta)
        meta->timestamp = rxmeta.timestamp.GetTicks();
    return samplesRead;
}

uint32_t TRXLooper::Receive(lime::complex32f_t* const* samples, uint32_t count, StreamRxMeta* meta)
{
    return StreamRxTemplate<complex32f_t>(samples, count, meta, chrono::microseconds(1000000));
}

uint32_t TRXLooper::Receive(lime::complex16_t* const* samples, uint32_t count, StreamRxMeta* meta)
{
    return StreamRxTemplate<complex16_t>(samples, count, meta, chrono::microseconds(1000000));
}

uint32_t TRXLooper::Receive(lime::complex12_t* const* samples, uint32_t count, StreamRxMeta* meta)
{
    return StreamRxTemplate<complex12_t>(samples, count, meta, chrono::microseconds(1000000));
}

OpStatus TRXLooper::TxSetup()
{
    OpStatus status = mTxArgs.dma->Initialize();
    if (status != OpStatus::Success)
        return status;

    mTx.samplesInPkt = defaultSamplesInPkt;
    mTx.terminate.store(false, std::memory_order_relaxed);

    mTx.lastTimestamp.store(0, std::memory_order_relaxed);
    const int chCount = std::max(mConfig.channels.at(lime::TRXDir::Rx).size(), mConfig.channels.at(lime::TRXDir::Tx).size());
    assert(chCount > 0);
    const int sampleSize = (mConfig.linkFormat == DataFormat::I16 ? 4 : 3); // sizeof IQ pair

    const GatewareFeatures gw = fpga->GetFeatures();
    uint32_t packetSize;
    if (gw.hasConfigurableStreamPacketSize)
    {
        mTx.samplesInPkt = 256 / chCount;
        if (mConfig.extraConfig.tx.samplesInPacket != 0)
        {
            mTx.samplesInPkt = mConfig.extraConfig.tx.samplesInPacket;
            lime::debug("Tx samples override %i", mTx.samplesInPkt);
        }
        packetSize = GetPacketSizeForBusSize(mTx.samplesInPkt, sampleSize, chCount, 32);
        const int headerSize = sizeof(StreamHeader);
        mTx.samplesInPkt = (packetSize - headerSize) / (sampleSize * chCount);
    }
    else
    {
        // FT601 USB encounters random BUS and IOMMU errors if transmitting not in 4096 byte chunks
        mTx.samplesInPkt = 4080 / sampleSize / chCount;
        packetSize = 4096;
    }

    mTx.packetsToBatch = 16; // Tx packets can be flushed early without filling whole batch
    // aim batch size to desired data output period, ~100us should be good enough
    if (mConfig.hintSampleRate > 0)
        mTx.packetsToBatch = std::floor((0.0005 * mConfig.hintSampleRate) / mTx.samplesInPkt);

    if (mConfig.extraConfig.tx.packetsInBatch != 0)
    {
        mTx.packetsToBatch = mConfig.extraConfig.tx.packetsInBatch;
    }

    const auto dmaChunks{ mTxArgs.dma->GetBuffers() };
    const auto dmaBufferSize = dmaChunks.front().size;

    mTx.packetsToBatch = std::clamp<uint8_t>(mTx.packetsToBatch, 1, dmaBufferSize / packetSize);

    std::vector<uint8_t*> dmaBuffers(dmaChunks.size());
    for (uint32_t i = 0; i < dmaChunks.size(); ++i)
    {
        dmaBuffers[i] = dmaChunks[i].buffer;
    }

    mTxArgs.buffers = std::move(dmaBuffers);
    mTxArgs.bufferSize = dmaBufferSize;
    mTxArgs.packetSize = packetSize;
    mTxArgs.packetsToBatch = mTx.packetsToBatch;
    mTxArgs.samplesInPacket = mTx.samplesInPkt;

    {
        float bufferTimeDuration;
        if (mConfig.hintSampleRate)
            bufferTimeDuration = mTx.samplesInPkt * mTx.packetsToBatch / mConfig.hintSampleRate;
        else
            bufferTimeDuration = 0;
        char msg[256];
        std::snprintf(msg,
            sizeof(msg),
            "Tx%i Setup: samplesInTxPkt:%i maxTxPktInBatch:%i, batchSizeInTime:%gus",
            chipId,
            mTx.samplesInPkt,
            mTx.packetsToBatch,
            bufferTimeDuration * 1e6);
        if (showStats)
            printf("%s\n", msg);
        if (mCallback_logMessage)
            mCallback_logMessage(LogLevel::Verbose, msg);
    }

    const int packetsInFIFO = 0.25 * mConfig.hintSampleRate / (mTx.packetsToBatch * mTx.samplesInPkt); // buffer 0.25 second of data
    mTx.packetsPool = std::make_unique<PacketsFIFO<StreamPacket*>>(packetsInFIFO);
    const uint32_t userSampleSize = mConfig.format == DataFormat::F32 ? sizeof(lime::complex32f_t) : sizeof(lime::complex16_t);
    for (uint32_t i = 0; i < mTx.packetsPool->max_size(); ++i)
        mTx.packetsPool->push(new StreamPacket(mTx.packetsToBatch * mTx.samplesInPkt, chCount, userSampleSize));
    mTx.fifo = std::make_unique<PacketsFIFO<StreamPacket*>>(packetsInFIFO);

    mTx.terminate.store(false, std::memory_order_relaxed);
    mTx.terminateWorker.store(false, std::memory_order_relaxed);
    auto TxLoopFunction = std::bind(&TRXLooper::TxWorkLoop, this);

    const auto schedulingPolicy = ThreadPolicy::REALTIME;
    mTx.thread = std::thread(TxLoopFunction);
    SetOSThreadPriority(ThreadPriority::HIGHEST, schedulingPolicy, &mTx.thread);
#ifdef __linux__
    char threadName[16]; // limited to 16 chars, including null byte.
    snprintf(threadName, sizeof(threadName), "lime:Tx%i", chipId);
    pthread_setname_np(mTx.thread.native_handle(), threadName);
#endif

    // Initialize DMA
    mTxArgs.dma->Enable(true);

    lime::debug("TxSetup wait for Tx worker.");
    // wait for Tx thread to be ready
    {
        std::unique_lock lck{ mTx.mutex };
        while (mTx.stage.load(std::memory_order_relaxed) < Stream::ReadyStage::WorkerReady)
            mTx.cv.wait(lck);
    }
    return OpStatus::Success;
}

void TRXLooper::TxWorkLoop()
{
    lime::debug("Tx worker thread ready.");
    // signal that thread is ready for work
    {
        std::unique_lock lck{ mTx.mutex };
        mTx.stage.store(Stream::ReadyStage::WorkerReady, std::memory_order_relaxed);
        mTx.cv.notify_all();
    }

    while (!mTx.terminateWorker.load(std::memory_order_relaxed))
    {
        // thread ready for work, just wait for stream enable
        {
            std::unique_lock lk{ streamMutex };
            while (!mStreamEnabled && !mTx.terminateWorker.load(std::memory_order_relaxed))
                streamActive.wait_for(lk, std::chrono::milliseconds(100));
        }
        if (!mStreamEnabled)
            continue;

        mTx.stage.store(Stream::ReadyStage::Active, std::memory_order_relaxed);
        TransmitPacketsLoop();
        {
            std::unique_lock lk{ mTx.mutex };
            mTx.stage.store(Stream::ReadyStage::WorkerReady, std::memory_order_relaxed);
            mTx.cv.notify_all();
        }
    }
    mTx.stage.store(Stream::ReadyStage::Disabled, std::memory_order_relaxed);
    lime::debug("Tx worker thread shutdown.");
}

static void TxPacketPadding(FPGA_TxDataPacket& packet, DataFormat linkFormat, uint8_t channelCount)
{
    // Tx data transfers have to be multiple of the bus size
    constexpr uint16_t busWidthBytes = 32;
    const int frameSize = (linkFormat == DataFormat::I12 ? 3 : 4) * channelCount;
    const int headerSize = sizeof(StreamHeader);
    const int maxPacketSize = 4096;

    uint16_t payloadSize = packet.GetPayloadSize();
    uint32_t samplesCount = payloadSize / frameSize;
    uint32_t packetSize = headerSize + payloadSize;

    while (packetSize % busWidthBytes)
    {
        packetSize += frameSize;
        ++samplesCount;

        if (packetSize >= maxPacketSize)
            break;
    }
    if (packetSize % busWidthBytes)
    {
        printf("Failed to pad Tx packet\n");
    }

    uint16_t paddingSize = packetSize - headerSize - packet.GetPayloadSize();
    if (paddingSize > 0)
    {
        std::memset(&packet.data[payloadSize], 0, paddingSize); // pad with zeroes
        packet.SetPayloadSize(payloadSize + paddingSize);
    }
}

void TRXLooper::TransmitPacketsLoop()
{
    lime::debug("Tx transmit loop start.");
    const bool isRxActive =
        true; // rx is always activated to provide timestamps // mConfig.channels.at(lime::TRXDir::Rx).size() > 0;
    const bool mimo = std::max(mConfig.channels.at(lime::TRXDir::Tx).size(), mConfig.channels.at(lime::TRXDir::Rx).size()) > 1;
    const bool compressed = mConfig.linkFormat == DataFormat::I12;
    constexpr int irqPeriod{ 4 }; // Interrupt request period

    const uint32_t bufferCount = mTxArgs.buffers.size();
    const std::vector<uint8_t*>& dmaBuffers{ mTxArgs.buffers };

    StreamStats& stats = mTx.stats;

    auto& fifo = mTx.fifo;

    int64_t totalBytesSent = 0; //for data rate calculation
    Timespec lastTS;

    struct PendingWrite {
        uint32_t id;
        uint8_t* data;
        uint32_t size;
    };
    std::queue<PendingWrite> pendingWrites;

    uint32_t stagingBufferIndex = 0;
    StreamPacket* srcPkt = nullptr;

    mTxArgs.dma->BufferOwnership(0, DataTransferDirection::DeviceToHost);

    bool outputReady = false;

    AvgRmsCounter txTSAdvance;

    auto t1{ std::chrono::steady_clock::now() };
    auto t2 = t1;

    DeltaVariable<int32_t> underrun(0);
    DeltaVariable<int32_t> loss(0);

    uint64_t lastHwIndex = 0;
    DMATransactionCounter counters;

    FPGA_TxDataPacket tempPacket;

    uint8_t* outputTail = dmaBuffers[0];
    uint64_t packetsCounter = 0;

    if (mConfig.timestampType == TimestampType::UNIX_EPOCH)
    {
        std::unique_lock lk{ startTimeMutex };
        while (!startUnixTimeSet && !mTx.terminateWorker.load(std::memory_order_relaxed))
        {
            startTimeIsSet.wait_for(lk, std::chrono::milliseconds(100));
        }
    }

    while (mTx.terminate.load(std::memory_order_relaxed) == false)
    {
        IDMA::State dma{ mTxArgs.dma->GetCounters() };
        int64_t counterDiff = ReadySlots(dma.transfersCompleted, lastHwIndex, 65536);
        lastHwIndex = dma.transfersCompleted;
        counters.completed += counterDiff;

        // process pending transactions
        while (!pendingWrites.empty() && counterDiff > 0)
        {
            PendingWrite& dataBlock = pendingWrites.front();
            totalBytesSent += dataBlock.size;
            stats.bytesTransferred += dataBlock.size;
            pendingWrites.pop();
            --counterDiff;
        }

        t2 = std::chrono::steady_clock::now();
        const auto timePeriod{ std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() };
        if (timePeriod >= statsPeriod_ms || mTx.terminate.load(std::memory_order_relaxed))
        {
            t1 = t2;
            double dataRate = 1000.0 * totalBytesSent / timePeriod;
            mTx.stats.dataRate_Bps = dataRate;

            double avgTxAdvance = 0, rmsTxAdvance = 0;
            txTSAdvance.GetResult(avgTxAdvance, rmsTxAdvance);
            loss.set(stats.loss);
            if (showStats || mCallback_logMessage)
            {
                char msg[512];
                std::snprintf(msg,
                    sizeof(msg) - 1,
                    "%s Tx%i: %3.3f MB/s | TS:%s pkt:%" PRIi64 " u:%i(%+i) l:%i(%+i) dma:%" PRIu64 "/%" PRIu64 "(%+" PRIi64
                    ") tsAdvance:%+.0f/%+.0f/%+.0f%s, f:%" PRIuPTR,
                    mTxArgs.dma->GetName().c_str(),
                    chipId,
                    dataRate / 1000000.0,
                    TimestampToString(lastTS, mConfig.timestampType).c_str(),
                    stats.packets,
                    underrun.value(),
                    underrun.delta(),
                    loss.value(),
                    loss.delta(),
                    counters.completed,
                    counters.requests,
                    counters.requests - counters.completed,
                    txTSAdvance.Min(),
                    avgTxAdvance,
                    txTSAdvance.Max(),
                    (mConfig.hintSampleRate ? "us" : ""),
                    fifo->size());
                if (showStats)
                    lime::info("%s", msg);
                if (mCallback_logMessage)
                {
                    bool showAsWarning = underrun.delta() || loss.delta();
                    LogLevel level = showAsWarning ? LogLevel::Warning : LogLevel::Debug;
                    mCallback_logMessage(level, msg);
                }
            }
            loss.checkpoint();
            underrun.checkpoint();
            totalBytesSent = 0;
        }

        bool reportProblems = false;
        // collect and transform samples data to output buffer
        while (!outputReady && !mTx.terminate.load(std::memory_order_relaxed))
        {
            if (!srcPkt)
            {
                if (!fifo->pop(&srcPkt, true, chrono::microseconds(100000)))
                {
                    std::this_thread::yield();
                    break;
                }
                if (mConfig.extraConfig.negateQ)
                    NegateQChannel(srcPkt, mConfig.format);

                lastTS = srcPkt->meta.timestamp;
                if (mConfig.timestampType == TimestampType::UNIX_EPOCH)
                {
                    srcPkt->meta.timestamp = srcPkt->meta.timestamp - Timespec(startUnixTime);
                    if (srcPkt->meta.useTimestamp && srcPkt->meta.timestamp.GetSeconds() < 0) // Drop packets that are in the past
                    {
                        reportProblems = true;
                        ++stats.underrun;
                        srcPkt->Reset();
                        mTx.packetsPool->push(srcPkt, true);
                        srcPkt = nullptr;
                        break;
                    }
                }
            }

            uint32_t payloadOffset = tempPacket.GetPayloadSize();
            uint8_t* payload = &tempPacket.data[payloadOffset];

            if (payloadOffset == 0)
            {
                tempPacket.ignoreTimestamp(!srcPkt->meta.useTimestamp);
                if (mConfig.timestampType == TimestampType::SAMPLE_TICKS)
                    tempPacket.counter = srcPkt->meta.timestamp.GetTicks();
                else
                {
                    tempPacket.counter = srcPkt->meta.timestamp.GetSeconds() << 32;
                    tempPacket.counter |=
                        uint64_t(srcPkt->meta.timestamp.GetFracSeconds() * mConfig.hintSampleRate * ticksPerSample);
                }
            }

            uint32_t bytesForFrame = (compressed ? 3 : 4) * (mimo ? 2 : 1);

            uint32_t samplesFilled = payloadOffset / bytesForFrame;
            uint32_t samplesToConsume = std::min(mTxArgs.samplesInPacket - samplesFilled, srcPkt->samples.size());

            DataConversion conversion;
            conversion.srcFormat = mConfig.format; //DataFormat::F32;
            conversion.destFormat = compressed ? DataFormat::I12 : DataFormat::I16;
            conversion.channelCount = mimo ? 2 : 1;

            int samplesDataSize = Interleave(payload, srcPkt->samples.front(), samplesToConsume, conversion);
            srcPkt->samples.pop(samplesToConsume);
            srcPkt->meta.timestamp.AddTicks(samplesToConsume * ticksPerSample);
            payloadOffset += samplesDataSize;
            tempPacket.SetPayloadSize(payloadOffset);
            assert(payloadOffset > 0);

            samplesFilled = payloadOffset / bytesForFrame;
            bool isPacketFull = samplesFilled == mTxArgs.samplesInPacket;
            bool doFlush = srcPkt->meta.flush && srcPkt->samples.size() == 0;
            if (isPacketFull || doFlush)
            {
                ++packetsCounter;

                TxPacketPadding(tempPacket, mConfig.linkFormat, conversion.channelCount);

                const int producedDataSize = sizeof(StreamHeader) + tempPacket.GetPayloadSize();
                memcpy(outputTail, &tempPacket, producedDataSize);
                outputTail += producedDataSize;
                tempPacket.ClearHeader();
            }

            doFlush |= packetsCounter == mTxArgs.packetsToBatch;

            if (srcPkt->samples.empty())
            {
                mTx.packetsPool->push(srcPkt, true);
                srcPkt = nullptr;
            }

            if (doFlush)
            {
                stats.packets += packetsCounter;
                packetsCounter = 0;
                outputReady = true;
                break;
            }
        }

        // one callback for the entire batch
        if (reportProblems && mConfig.statusCallback)
            mConfig.statusCallback(true, &stats, mConfig.userData);

        bool canSend = pendingWrites.size() < bufferCount - 1;
        if (!canSend)
        {
            if (mConfig.extraConfig.usePoll)
                canSend = mTxArgs.dma->Wait() == OpStatus::Success;
            else
                std::this_thread::yield();
            continue;
        }

        if (!outputReady)
            continue;

        FPGA_TxDataPacket* pkt = reinterpret_cast<FPGA_TxDataPacket*>(dmaBuffers[stagingBufferIndex]);
        if (!pkt->getIgnoreTimestamp() && isRxActive) // Rx is needed for current timestamp
        {
            int64_t rxNow = mRx.lastTimestamp.load(std::memory_order_relaxed);
            const int64_t txAdvance = pkt->counter - rxNow;
            if (mConfig.hintSampleRate)
            {
                int64_t timeAdvance = ts_to_us(mConfig.hintSampleRate, txAdvance);
                txTSAdvance.Add(timeAdvance);
            }
            else
                txTSAdvance.Add(txAdvance);
            if (txAdvance <= 0)
            {
                reportProblems = true;
                underrun.add(1);
                ++stats.underrun;
            }
        }

        if (reportProblems && mConfig.statusCallback)
            mConfig.statusCallback(true, &stats, mConfig.userData);

        uint32_t bytesToSend = outputTail - dmaBuffers[stagingBufferIndex];
        PendingWrite wrInfo{ stagingBufferIndex, dmaBuffers[stagingBufferIndex], bytesToSend };
        bool requestIRQ = (wrInfo.id % irqPeriod) == 0;
        // DMA memory is write only, to read from the buffer will trigger Bus errors
        mTxArgs.dma->BufferOwnership(stagingBufferIndex, DataTransferDirection::HostToDevice);
        const OpStatus status{ mTxArgs.dma->SubmitRequest(
            stagingBufferIndex, wrInfo.size, DataTransferDirection::HostToDevice, requestIRQ) };
        if (status != OpStatus::Success)
        {
            lime::error("Failed to submit dma write");
            ++stats.overrun;
            mTxArgs.dma->Wait();
            continue;
        }

        pendingWrites.push(wrInfo);
        stagingBufferIndex = (stagingBufferIndex + 1) % bufferCount;
        mTxArgs.dma->BufferOwnership(stagingBufferIndex, DataTransferDirection::DeviceToHost);
        ++counters.requests;

        outputReady = false;
        stats.timestamp = lastTS.GetTicks();
        outputTail = dmaBuffers[stagingBufferIndex];
    }
    lime::debug("Tx transmit loop end.");
}

void TRXLooper::TxTeardown()
{
    if (mTx.stage.load(std::memory_order_relaxed) != Stream::ReadyStage::Disabled)
    {
        lime::debug("TxTeardown wait for Tx worker shutdown.");
        mTx.terminateWorker.store(true, std::memory_order_relaxed);
        {
            std::unique_lock lck{ streamMutex };
            mTx.terminateWorker.store(true, std::memory_order_relaxed);
            streamActive.notify_all();
        }

        mTx.terminate.store(true, std::memory_order_relaxed);
        try
        {
            if (mTx.thread.joinable())
                mTx.thread.join();
        } catch (...)
        {
            lime::error("Failed to join TRXLooper Tx thread"s);
        }
    }

    if (mTx.stagingPacket)
    {
        delete mTx.stagingPacket;
        mTx.stagingPacket = nullptr;
    }
    if (mTx.packetsPool)
    {
        while (mTx.packetsPool->pop(&mTx.stagingPacket, false))
            delete mTx.stagingPacket;
        mTx.stagingPacket = nullptr;
    }

    delete mTx.fifo.release();
    delete mTx.packetsPool.release();
}

template<class T>
uint32_t TRXLooper::StreamMetaToStreamTxMeta(
    const T* const* samples, uint32_t count, const StreamMeta* meta, std::chrono::microseconds timeout)
{
    StreamTxMeta txmeta;
    txmeta.hasTimestamp = meta ? meta->waitForTimestamp : false;
    txmeta.flags = meta ? (meta->flushPartialPacket ? StreamTxMeta::EndOfBurst : 0) : 0;
    if (txmeta.hasTimestamp)
    {
        if (mConfig.timestampType == TimestampType::SAMPLE_TICKS)
            txmeta.timestamp = Timespec(meta->timestamp / mConfig.hintSampleRate);
        else
            txmeta.timestamp = Timespec(meta->timestamp >> 32, (meta->timestamp & 0xFFFFFFFF) / 1e9);
    }
    return StreamTxTemplate(samples, count, &txmeta, timeout);
}

template<class T>
uint32_t TRXLooper::StreamTxTemplate(
    const T* const* samples, uint32_t count, const StreamTxMeta* meta, chrono::microseconds timeout)
{
    const bool useChannelB = mConfig.channels.at(lime::TRXDir::Tx).size() > 1;
    const bool useTimestamp = meta ? (meta->hasTimestamp) : false;
    const bool flush = meta ? (meta->flags & StreamTxMeta::EndOfBurst) : false;

    Timespec ts = meta->timestamp;
    ts.SetTickRate(ticksPerSample * mConfig.hintSampleRate);

    uint32_t samplesRemaining = count;

    bool timeGap = true; // expectedTS != lime::Timespec(meta->timestamp);

    if (mTx.stagingPacket && timeGap)
    {
        if (!mTx.fifo->push(mTx.stagingPacket, true, timeout))
            return 0;

        mTx.stagingPacket = nullptr;
    }

    assert(samples);
    assert(samples[0]);
    if (useChannelB)
        assert(samples[1]);
    const T* src[2] = { samples[0], useChannelB ? samples[1] : nullptr };
    while (samplesRemaining > 0)
    {
        if (!mTx.stagingPacket)
        {
            mTx.packetsPool->pop(&mTx.stagingPacket, true);
            if (!mTx.stagingPacket)
                break;

            mTx.stagingPacket->Reset();
            mTx.stagingPacket->meta.timestamp = ts;
            mTx.stagingPacket->meta.useTimestamp = useTimestamp;
        }

        int consumed = mTx.stagingPacket->samples.push(src, samplesRemaining);
        src[0] += consumed;
        if (useChannelB)
            src[1] += consumed;

        samplesRemaining -= consumed;
        ts.AddTicks(consumed * ticksPerSample);

        bool pushPacket = mTx.stagingPacket->samples.isFull();

        if (samplesRemaining == 0 && flush)
        {
            mTx.stagingPacket->meta.flush = flush;
            pushPacket = true;
        }

        if (pushPacket)
        {
            if (!mTx.fifo->push(mTx.stagingPacket, true, chrono::microseconds(1000000)))
                break;

            mTx.stagingPacket = nullptr;
        }
    }

    return count - samplesRemaining;
}

/// @brief Transmits packets from from this specific stream.
/// @param samples The buffer of the samples to transmit.
/// @param count The amount of samples to transmit.
/// @param meta The metadata of the packets of the stream.
/// @return The amount of samples transmitted.
uint32_t TRXLooper::StreamTx(
    const lime::complex32f_t* const* samples, uint32_t count, const StreamMeta* meta, std::chrono::microseconds timeout)
{
    return StreamMetaToStreamTxMeta(samples, count, meta, timeout);
}

/// @copydoc TRXLooper::StreamTx()
uint32_t TRXLooper::StreamTx(
    const lime::complex16_t* const* samples, uint32_t count, const StreamMeta* meta, std::chrono::microseconds timeout)
{
    return StreamMetaToStreamTxMeta(samples, count, meta, timeout);
}

/// @copydoc TRXLooper::StreamTx()
uint32_t TRXLooper::StreamTx(
    const lime::complex12_t* const* samples, uint32_t count, const StreamMeta* meta, std::chrono::microseconds timeout)
{
    return StreamMetaToStreamTxMeta(samples, count, meta, timeout);
}

uint32_t TRXLooper::Transmit(const lime::complex32f_t* const* samples, uint32_t count, const StreamTxMeta* meta)
{
    return StreamTxTemplate(samples, count, meta, chrono::microseconds(100000));
}

uint32_t TRXLooper::Transmit(const lime::complex16_t* const* samples, uint32_t count, const StreamTxMeta* meta)
{
    return StreamTxTemplate(samples, count, meta, chrono::microseconds(100000));
}

uint32_t TRXLooper::Transmit(const lime::complex12_t* const* samples, uint32_t count, const StreamTxMeta* meta)
{
    return StreamTxTemplate(samples, count, meta, chrono::microseconds(100000));
}

/// @brief Gets Rx/Tx data transfer statistics.
/// @param rxStats Pointer to rx statistics structure, (Optional, can be NULL)
/// @param txStats Pointer to rx statistics structure, (Optional, can be NULL)
void TRXLooper::StreamStatus(StreamStats* rxStats, StreamStats* txStats)
{
    if (txStats)
    {
        *txStats = mTx.stats;
        if (mTx.fifo)
            txStats->FIFO = { mTx.fifo->max_size(), mTx.fifo->size() };
        else
            txStats->FIFO = { 1, 0 };
    }
    if (rxStats)
    {
        *rxStats = mRx.stats;
        if (mRx.fifo)
            rxStats->FIFO = { mRx.fifo->max_size(), mRx.fifo->size() };
        else
            rxStats->FIFO = { 1, 0 };
    }
}

/// @copydoc SDRDevice::UploadTxWaveform()
/// @param fpga The FPGA device to use.
/// @param dma The data channel to use.
OpStatus TRXLooper::UploadTxWaveform(FPGA* fpga,
    std::shared_ptr<IDMA> dma,
    const lime::StreamConfig& config,
    uint8_t moduleIndex,
    const void** samples,
    uint32_t count)
{
    const int samplesInPkt = 256;
    const bool useChannelB = config.channels.at(lime::TRXDir::Tx).size() > 1;
    const bool mimo = config.channels.at(lime::TRXDir::Tx).size() == 2;
    const bool compressed = config.linkFormat == DataFormat::I12;
    fpga->SelectModule(moduleIndex);
    fpga->WriteRegister(0x000C, mimo ? 0x3 : 0x1); //channels 0,1
    if (config.linkFormat == DataFormat::I16)
        fpga->WriteRegister(0x000E, 0x0); //16bit samples
    else
        fpga->WriteRegister(0x000E, 0x2); //12bit samples

    fpga->WriteRegister(0x000D, 0x4); // WFM_LOAD

    DataConversion conversion;
    conversion.srcFormat = config.format;
    conversion.destFormat = compressed ? DataFormat::I12 : DataFormat::I16;
    conversion.channelCount = mimo ? 2 : 1;

    const auto dmaChunks{ dma->GetBuffers() };

    std::vector<uint8_t*> dmaBuffers(dmaChunks.size());
    for (uint32_t i = 0; i < dmaChunks.size(); ++i)
        dmaBuffers[i] = dmaChunks[i].buffer;

    dma->Enable(true);

    uint32_t samplesRemaining = count;
    uint8_t stagingBufferIndex = 0;

    const uint8_t* src[2] = { static_cast<const uint8_t*>(samples[0]),
        static_cast<const uint8_t*>(useChannelB ? samples[1] : nullptr) };
    const int sampleSize = config.format == DataFormat::F32 ? sizeof(lime::complex32f_t) : sizeof(lime::complex16_t);

    while (samplesRemaining > 0)
    {
        // IDMA::State state{ dma->GetCounters() };
        dma->Wait(); // block until there is a free DMA buffer

        int samplesToSend = samplesRemaining > samplesInPkt ? samplesInPkt : samplesRemaining;
        int samplesDataSize = 0;

        dma->BufferOwnership(stagingBufferIndex, DataTransferDirection::DeviceToHost);

        FPGA_TxDataPacket* pkt = reinterpret_cast<FPGA_TxDataPacket*>(dmaBuffers[stagingBufferIndex]);
        pkt->counter = 0;
        pkt->reserved[0] = 0;

        samplesDataSize = Interleave(pkt->data, reinterpret_cast<const void**>(src), samplesToSend, conversion);

        src[0] += samplesToSend * sampleSize;
        if (useChannelB)
            src[1] += samplesToSend * sampleSize;

        int payloadSize = (samplesDataSize / 4) * 4;
        if (samplesDataSize % 4 != 0)
            lime::warning("Packet samples count not multiple of 4"s);
        pkt->reserved[2] = (payloadSize >> 8) & 0xFF; //WFM loading
        pkt->reserved[1] = payloadSize & 0xFF; //WFM loading
        pkt->reserved[0] = 0x1 << 5; //WFM loading

        dma->BufferOwnership(stagingBufferIndex, DataTransferDirection::HostToDevice);

        size_t transferSize = 16 + payloadSize;
        assert(transferSize <= dmaChunks.front().size);
        bool requestIRQ = (stagingBufferIndex % 4) == 0;

        // DMA memory is write only, to read from the buffer will trigger Bus errors
        const OpStatus status{ dma->SubmitRequest(
            stagingBufferIndex, transferSize, DataTransferDirection::HostToDevice, requestIRQ) };
        if (status != OpStatus::Success)
        {
            dma->Enable(false);
            return ReportError(OpStatus::IOFailure, "Failed to submit dma write (%i) %s", errno, strerror(errno));
        }

        samplesRemaining -= samplesToSend;
        stagingBufferIndex = (stagingBufferIndex + 1) % dmaBuffers.size();
    }

    // Give some time to load samples to FPGA
    std::this_thread::sleep_for(std::chrono::milliseconds(500));
    dma->Enable(false);

    fpga->StopWaveformPlayback();
    if (samplesRemaining != 0)
        return ReportError(OpStatus::Error, "Failed to upload waveform"s);

    return OpStatus::Success;
}

} // namespace lime
