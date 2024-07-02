import struct
import pandas as pd

from pick import pick


class EDFHandler:
    def __init__(self, edf_file_path):
        self.edf_file_path = edf_file_path

        self.edf_file = self.load_edf_file()
        self.header_info = self.get_header_info()

    def load_edf_file(self):
        return open(self.edf_file_path, "rb").read()

    def fix_edf_file(self):
        chars_replacements = {
            b"\xDF": "s",
            b"\xE4": "a",
            b"\xC4": "A",
            b"\xF6": "o",
            b"\xD6": "O",
            b"\xFC": "u",
            b"\xDC": "U",
            b"\xb0": "d",  # degrees symbol
        }

        for byte_val, replacement in chars_replacements.items():
            self.edf_file = self.edf_file.replace(
                byte_val, replacement.encode("latin-1")
            )

    def get_header_info(self):
        header = self.edf_file[:256].decode("latin-1")
        header_dict = {
            "version": header[:8].strip(),
            "patient_id": header[8:88].strip(),
            "recording_id": header[88:168].strip(),
            "start_date": header[168:176].strip(),
            "start_time": header[176:184].strip(),
            "header_bytes": header[184:192].strip(),
            "num_records": header[236:244].strip(),
            "duration": header[244:252].strip(),
            "num_signals": int(header[252:256].strip()),
        }

        return header_dict

    def get_signals_info(self):
        signals_bytes_chunk = self.edf_file[
            256 : 256 + self.header_info["num_signals"] * 16
        ]
        transducers_bytes_chunk = self.edf_file[
            256
            + (16 * self.header_info["num_signals"]) : 256
            + (96 * self.header_info["num_signals"])
        ]
        units_bytes_chunk = self.edf_file[
            256
            + (96 * self.header_info["num_signals"]) : 256
            + (104 * self.header_info["num_signals"])
        ]
        physical_min_bytes_chunk = self.edf_file[
            256
            + (104 * self.header_info["num_signals"]) : 256
            + (112 * self.header_info["num_signals"])
        ]
        physical_max_bytes_chunk = self.edf_file[
            256
            + (112 * self.header_info["num_signals"]) : 256
            + (120 * self.header_info["num_signals"])
        ]
        digital_min_bytes_chunk = self.edf_file[
            256
            + (120 * self.header_info["num_signals"]) : 256
            + (128 * self.header_info["num_signals"])
        ]
        digital_max_bytes_chunk = self.edf_file[
            256
            + (128 * self.header_info["num_signals"]) : 256
            + (136 * self.header_info["num_signals"])
        ]
        prefiltering_bytes_chunk = self.edf_file[
            256
            + (136 * self.header_info["num_signals"]) : 256
            + (216 * self.header_info["num_signals"])
        ]
        num_samples_bytes_chunk = self.edf_file[
            256
            + (216 * self.header_info["num_signals"]) : 256
            + (224 * self.header_info["num_signals"])
        ]
        reserved_bytes_chunk = self.edf_file[
            256
            + (224 * self.header_info["num_signals"]) : 256
            + (256 * self.header_info["num_signals"])
        ]

        signals_labels = [
            signals_bytes_chunk[i * 16 : i * 16 + 16].decode("latin-1").strip()
            for i in range(self.header_info["num_signals"])
        ]
        transducers_labels = [
            transducers_bytes_chunk[i * 80 : i * 80 + 80].decode("latin-1").strip()
            for i in range(self.header_info["num_signals"])
        ]
        units_labels = [
            units_bytes_chunk[i * 8 : i * 8 + 8].decode("latin-1").strip()
            for i in range(self.header_info["num_signals"])
        ]
        physical_min_values = [
            float(physical_min_bytes_chunk[i * 8 : i * 8 + 8].decode("latin-1").strip())
            for i in range(self.header_info["num_signals"])
        ]
        physical_max_values = [
            float(physical_max_bytes_chunk[i * 8 : i * 8 + 8].decode("latin-1").strip())
            for i in range(self.header_info["num_signals"])
        ]
        digital_min_values = [
            int(digital_min_bytes_chunk[i * 8 : i * 8 + 8].decode("latin-1").strip())
            for i in range(self.header_info["num_signals"])
        ]
        digital_max_values = [
            int(digital_max_bytes_chunk[i * 8 : i * 8 + 8].decode("latin-1").strip())
            for i in range(self.header_info["num_signals"])
        ]
        prefiltering_values = [
            prefiltering_bytes_chunk[i * 80 : i * 80 + 80].decode("latin-1").strip()
            for i in range(self.header_info["num_signals"])
        ]
        num_samples_per_record = [
            int(num_samples_bytes_chunk[i * 8 : i * 8 + 8].decode("latin-1").strip())
            for i in range(self.header_info["num_signals"])
        ]
        reserved_values = [
            reserved_bytes_chunk[i * 32 : i * 32 + 32].decode("latin-1").strip()
            for i in range(self.header_info["num_signals"])
        ]

        signals_header_info = {
            "signals_labels": signals_labels,
            "transducers_labels": transducers_labels,
            "units_labels": units_labels,
            "physical_min_values": physical_min_values,
            "physical_max_values": physical_max_values,
            "digital_min_values": digital_min_values,
            "digital_max_values": digital_max_values,
            "prefiltering_values": prefiltering_values,
            "num_samples_per_record": num_samples_per_record,
            "reserved_values": reserved_values,
        }

        return signals_header_info

    def read_record_chn(self, signal_name):
        header_info = self.get_header_info()
        signals_info = self.get_signals_info()

        n_records = int(header_info["num_records"])
        header_bytes = int(header_info["header_bytes"])

        signal_index = signals_info["signals_labels"].index(signal_name)
        n_samples_for_signal = signals_info["num_samples_per_record"][signal_index]

        records_data = []

        bytes_per_sample = 2
        total_samples_per_record = sum(signals_info["num_samples_per_record"])
        bytes_per_record = total_samples_per_record * bytes_per_sample

        offset = header_bytes
        for j in range(signal_index):
            offset += signals_info["num_samples_per_record"][j] * 2

        with open(self.edf_file_path, "rb") as file:
            file.seek(offset)

            for _ in range(n_records):
                samples = []
                for _ in range(n_samples_for_signal):
                    samples.append(struct.unpack("<h", file.read(2))[0])
                records_data.append(samples)

                # skip other signals
                bytes_to_skip = bytes_per_record - n_samples_for_signal * 2
                file.seek(bytes_to_skip, 1)

        flatten_records = [sample for record in records_data for sample in record]

        print(f"Read {signal_name} signal.")
        return flatten_records

    def export_signals_to_csv(self):
        label_names = self.get_signals_info()["signals_labels"]
        tranducer_names = self.get_signals_info()["transducers_labels"]
        unit_names = self.get_signals_info()["units_labels"]

        signal_names = [f'{label} - {tranducer} ({unit})' for label, tranducer, unit in zip(label_names, tranducer_names, unit_names)]
        signal_names.append("All signals")

        title = "Select signals to export to CSV: "
        selected_signals = pick(
            signal_names, title, multiselect=True, min_selection_count=1
        )
        selected_signals = [signal_name for signal_name, _ in selected_signals]
        selected_signals = [signal_name.split(" - ")[0] for signal_name in selected_signals]

        if "All signals" in selected_signals:
            selected_signals = signal_names[:-1]

        signals_data = []
        for signal_name in selected_signals:
            signal_data = self.read_record_chn(signal_name)
            signals_data.append(signal_data)

        df = pd.DataFrame(signals_data).T
        df.columns = selected_signals

        output_path = self.edf_file_path.replace(".edf", ".csv")
        df.to_csv(output_path, index=False)
        print(f"Signals exported to: {output_path}")

        return


edf_file = input("Enter the EDF file path: ")
edf = EDFHandler(edf_file)
edf.export_signals_to_csv()
