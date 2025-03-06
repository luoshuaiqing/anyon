#!/usr/bin/env python
# coding: utf-8

import matplotlib.pyplot as plt
import csv
import numpy as np
import glob
import os
import pandas as pd
import re
import math
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import LinearSegmentedColormap, Normalize
import time
import datetime
from scipy import signal
from scipy.signal import find_peaks
from scipy.linalg import sqrtm

RED = "\033[91m"
GREEN = "\033[92m"
RESET = "\033[0m"  # Reset to default color


def combine_binary(bin_list):
    combined_binary = "".join(bin_list)
    return combined_binary


def fortyeightbit_change_o_single(
    result1, ave_i_0, ave_q_0, color="blue", color1="red"
):
    result1 = np.array(result1)
    indices = np.where(result1 < 0)
    result1[indices] += 65536
    binary_array = [bin(decimal)[2:].zfill(16) for decimal in result1]
    combined_array_i = []
    combined_array_q = []
    for i in range(85):
        binary_parts_i = [
            binary_array[i * 3 + 1],
            binary_array[i * 3 + 2],
            binary_array[i * 3 + 3],
        ]
        binary_parts_q = [
            binary_array[i * 3 + 257],
            binary_array[i * 3 + 258],
            binary_array[i * 3 + 259],
        ]
        combined_array_i.append(combine_binary(binary_parts_i))
        combined_array_q.append(combine_binary(binary_parts_q))

    decimal_array_i = []
    for binary_str in combined_array_i:
        if binary_str[0] == "1":
            decimal_value = int(binary_str, 2) - 2**48
        else:
            decimal_value = int(binary_str, 2)
        decimal_array_i.append(decimal_value)
    decimal_array_q = []
    for binary_str in combined_array_q:
        if binary_str[0] == "1":
            decimal_value = int(binary_str, 2) - 2**48
        else:
            decimal_value = int(binary_str, 2)
        decimal_array_q.append(decimal_value)

    single_array_ampl = []
    single_array_phase = []
    for i in range(85):
        single_array_ampl.append(
            math.sqrt(decimal_array_i[i] ** 2 + decimal_array_q[i] ** 2)
        )
        single_array_phase.append(math.atan2(decimal_array_i[i], decimal_array_q[i]))

    return single_array_ampl, single_array_phase, decimal_array_i, decimal_array_q


def process_file(file_path, ave_i_pretest, ave_q_pretest):
    """Process a single file to calculate ampl and phase."""
    df = pd.read_csv(file_path, header=None, on_bad_lines="skip")
    data_arrays = [np.array(row) for _, row in df.iterrows()]
    scatter = np.array(data_arrays)
    single_file_result = fortyeightbit_change_o_single(
        scatter[0], ave_i_pretest, ave_q_pretest
    )
    ampl_85 = single_file_result[0]
    phase_85 = single_file_result[1]
    i_85 = single_file_result[2]
    q_85 = single_file_result[3]
    return ampl_85, phase_85, i_85, q_85


def process_folder(folder_path, ave_i_pretest, ave_q_pretest, keyword=""):
    """Process all files in a folder and calculate ampl and phase for each."""
    files = [
        os.path.join(folder_path, f)
        for f in os.listdir(folder_path)
        if f.endswith(".csv") and keyword in f
    ]
    ampl_list = []
    phase_list = []
    i_list = []
    q_list = []

    for file in files:
        ampl, phase, i, q = process_file(file, ave_i_pretest, ave_q_pretest)
        ampl_list.append(ampl)
        phase_list.append(phase)
        i_list.append(i)
        q_list.append(q)

    return ampl_list, phase_list, i_list, q_list


folder_path_q1 = "test2_C12bell_freq72.875_EQ2RQ1RQ2_campl_50_DC_0.2mA_Z/readout1/test"  # Replace with the actual folder path
folder_path_q2 = "test2_C12bell_freq72.875_EQ2RQ1RQ2_campl_50_DC_0.2mA_Z/readout2/test2"  # Replace with the actual folder path
ave_i_pretest_q1 = 0  #
ave_q_pretest_q1 = 0  #
ave_i_pretest_q2 = 0  #
ave_q_pretest_q2 = 0

ampl_list_q1, phase_list_q1, i_list_q1, q_list_q1 = process_folder(
    folder_path_q1, ave_i_pretest_q1, ave_q_pretest_q1
)
ampl_list_q2, phase_list_q2, i_list_q2, q_list_q2 = process_folder(
    folder_path_q2, ave_i_pretest_q2, ave_q_pretest_q2
)

plt.figure(figsize=(12, 6))

plt.subplot(1, 2, 1)
for i in range(len(i_list_q1)):
    plt.scatter(i_list_q1[i], q_list_q1[i], label=f"File {i+1}", alpha=0.6)
plt.xlabel("I Quadrature")
plt.ylabel("Q Quadrature")
plt.title("IQ Plane Scatter Plot - Qubit 1")
plt.legend()
plt.grid(True)

plt.subplot(1, 2, 2)
for i in range(len(i_list_q2)):
    plt.scatter(i_list_q2[i], q_list_q2[i], label=f"File {i+1}", alpha=0.6)
plt.xlabel("I Quadrature")
plt.ylabel("Q Quadrature")
plt.title("IQ Plane Scatter Plot - Qubit 2")
plt.legend()
plt.grid(True)

plt.tight_layout()
# plt.show()

# rotate to I
theta_rotate_q1 = 48
theta_rotate_q2 = 40


def rotate_iq(i_data, q_data, theta_deg):
    theta = np.radians(theta_deg)  # degree to rad
    i_data = np.array(i_data)  # Convert to numpy array
    q_data = np.array(q_data)  # Convert to numpy array
    i_rot = i_data * np.cos(theta) - q_data * np.sin(theta)
    q_rot = i_data * np.sin(theta) + q_data * np.cos(theta)
    return i_rot, q_rot


i_rotated_q1, q_rotated_q1 = rotate_iq(i_list_q1, q_list_q1, theta_rotate_q1)
i_rotated_q2, q_rotated_q2 = rotate_iq(i_list_q2, q_list_q2, theta_rotate_q2)

plt.figure(figsize=(12, 6))

plt.subplot(1, 2, 1)
for i in range(len(i_rotated_q1)):
    plt.scatter(i_rotated_q1[i], q_rotated_q1[i], label=f"File {i+1}", alpha=0.6)
plt.xlabel("I Quadrature")
plt.ylabel("Q Quadrature")
plt.title("IQ Plane Scatter Plot - Qubit 1 - Rotated")
plt.legend()
plt.grid(True)

plt.subplot(1, 2, 2)
for i in range(len(i_rotated_q2)):
    plt.scatter(i_rotated_q2[i], q_rotated_q2[i], label=f"File {i+1}", alpha=0.6)
plt.xlabel("I Quadrature")
plt.ylabel("Q Quadrature")
plt.title("IQ Plane Scatter Plot - Qubit 2 - Rotated")
plt.legend()
plt.grid(True)

plt.tight_layout()
# plt.show()


# set threshold
threshold_q1 = -0.665e9
threshold_q2 = 7.555e8

i_rotated_q1_1d = i_rotated_q1.flatten()
i_rotated_q2_1d = i_rotated_q2.flatten()


def analyze_qubit_states(i_rotated_q1, i_rotated_q2, threshold_q1, threshold_q2):
    counts = {"00": 0, "01": 0, "10": 0, "11": 0}

    for i1, i2 in zip(i_rotated_q1, i_rotated_q2):
        if i1 < threshold_q1 and i2 < threshold_q2:
            counts["00"] += 1
        elif i1 < threshold_q1 and i2 >= threshold_q2:
            counts["01"] += 1
        elif i1 >= threshold_q1 and i2 < threshold_q2:
            counts["10"] += 1
        else:
            counts["11"] += 1

    total = sum(counts.values())
    probabilities = {state: count / total for state, count in counts.items()}

    return probabilities


probabilities = analyze_qubit_states(
    i_rotated_q1_1d, i_rotated_q2_1d, threshold_q1, threshold_q2
)
print(probabilities)

plt.figure(figsize=(6, 4))
plt.bar(
    probabilities.keys(),
    probabilities.values(),
    color=["blue", "orange", "green", "red"],
)

plt.xlabel("Quantum States")
plt.ylabel("Probability")
plt.title("Quantum State Population Distribution")
plt.ylim(0, 1)
plt.grid(axis="y", linestyle="--", alpha=0.7)

for state, prob in probabilities.items():
    plt.text(state, prob + 0.02, f"{prob:.2f}", ha="center", fontsize=10)

# plt.show()


# 2q tomo
folder_path_q1 = "test2_C12bell_freq72.875_EQ2RQ1RQ2_campl_50_DC_0.2mA_tomo/readout1"  #
folder_path_q2 = "test2_C12bell_freq72.875_EQ2RQ1RQ2_campl_50_DC_0.2mA_tomo/readout2"  #

probs = np.zeros((4, 16))
for i in range(16):
    ampl_list_q1, phase_list_q1, i_list_q1, q_list_q1 = process_folder(
        folder_path_q1, ave_i_pretest_q1, ave_q_pretest_q1, f"_proj{i}"
    )
    ampl_list_q2, phase_list_q2, i_list_q2, q_list_q2 = process_folder(
        folder_path_q2, ave_i_pretest_q2, ave_q_pretest_q2, f"_proj{i}"
    )
    i_rotated_q1, q_rotated_q1 = rotate_iq(i_list_q1, q_list_q1, theta_rotate_q1)
    i_rotated_q2, q_rotated_q2 = rotate_iq(i_list_q2, q_list_q2, theta_rotate_q2)
    i_rotated_q1_1d = i_rotated_q1.flatten()
    i_rotated_q2_1d = i_rotated_q2.flatten()
    probabilities = analyze_qubit_states(
        i_rotated_q1_1d, i_rotated_q2_1d, threshold_q1, threshold_q2
    )
    probs[:, i] = [
        probabilities["00"],
        probabilities["01"],
        probabilities["10"],
        probabilities["11"],
    ]

stokes = np.zeros(15)
for i in range(1, 16):
    if i in [1, 2, 3]:
        stokes[i - 1] = (
            probs[0][i - 1] + probs[2][i - 1] - probs[1][i - 1] - probs[3][i - 1]
        )
    elif i in [4, 8, 12]:
        stokes[i - 1] = (
            probs[0][i - 1] + probs[1][i - 1] - probs[2][i - 1] - probs[3][i - 1]
        )
    elif i in [5, 6, 7, 9, 10, 11, 13, 14, 15]:
        stokes[i - 1] = (
            probs[0][i - 1] + probs[3][i - 1] - probs[1][i - 1] - probs[2][i - 1]
        )

# Derive the density matrix
I = np.array([[1, 0], [0, 1]])
sigma_x = np.array([[0, 1], [1, 0]])
sigma_y = np.array([[0, -1j], [1j, 0]])
sigma_z = np.array([[1, 0], [0, -1]])

# Density matrix Pauli operator basis
II = np.kron(I, I)
IX = np.kron(I, sigma_x)
IY = np.kron(I, sigma_y)
IZ = np.kron(I, sigma_z)
XI = np.kron(sigma_x, I)
XX = np.kron(sigma_x, sigma_x)
XY = np.kron(sigma_x, sigma_y)
XZ = np.kron(sigma_x, sigma_z)
YI = np.kron(sigma_y, I)
YX = np.kron(sigma_y, sigma_x)
YY = np.kron(sigma_y, sigma_y)
YZ = np.kron(sigma_y, sigma_z)
ZI = np.kron(sigma_z, I)
ZX = np.kron(sigma_z, sigma_x)
ZY = np.kron(sigma_z, sigma_y)
ZZ = np.kron(sigma_z, sigma_z)

rho = 0.25 * (
    II
    + stokes[0] * IX
    + stokes[1] * IY
    + stokes[2] * IZ
    + stokes[3] * XI
    + stokes[4] * XX
    + stokes[5] * XY
    + stokes[6] * XZ
    + stokes[7] * YI
    + stokes[8] * YX
    + stokes[9] * YY
    + stokes[10] * YZ
    + stokes[11] * ZI
    + stokes[12] * ZX
    + stokes[13] * ZY
    + stokes[14] * ZZ
)


print(f"The density matrix is:\n{np.round(rho, decimals=3)}")

ideal_state = np.array(
    [[0.25, 0, 0, 0], [0, 0.25, 0, 0], [0, 0, 0.25, 0], [0, 0, 0, 0.25]]
)
sqrt_ideal_state = sqrtm(ideal_state)
state_fidelity = (np.abs(sqrtm(sqrt_ideal_state @ rho @ sqrt_ideal_state).trace())) ** 2

print(f"The state fidelity is: {np.round(state_fidelity, decimals=4)}")


# 定义密度矩阵，基底顺序为 |00>, |01>, |10>, |11>
# rho = 0.5 * np.array([
#     [0,    0,   0,  0],
#     [0,    1, -1j, 0],
#     [0,  1j,   1,  0],
#     [0,    0,   0,  0]
# ])

rho_real = np.real(rho)
rho_imag = np.imag(rho)

plt.figure(figsize=(10, 4))

plt.subplot(1, 2, 1)
im1 = plt.imshow(rho_real, cmap="viridis", interpolation="none")
plt.title("Real Part")
plt.xticks(range(4), ["|00>", "|01>", "|10>", "|11>"])
plt.yticks(range(4), ["|00>", "|01>", "|10>", "|11>"])
plt.colorbar(im1, fraction=0.046, pad=0.04)

plt.subplot(1, 2, 2)
im2 = plt.imshow(rho_imag, cmap="viridis", interpolation="none")
plt.title("Imaginary Part")
plt.xticks(range(4), ["|00>", "|01>", "|10>", "|11>"])
plt.yticks(range(4), ["|00>", "|01>", "|10>", "|11>"])
plt.colorbar(im2, fraction=0.046, pad=0.04)

plt.tight_layout()
# plt.show()


def read_iq_from_csv(csv_path):
    df = pd.read_csv(csv_path)
    # Create a list to store all I and Q value arrays
    i_arrays = []
    q_arrays = []

    # For each of the 30 I/Q pairs
    for x in range(30):
        i_values = df.iloc[:, x].to_numpy()
        q_values = df.iloc[:, x + 30].to_numpy()
        i_arrays.append(i_values)
        q_arrays.append(q_values)

    return i_arrays, q_arrays


def generate_raw_data_np(i_values, q_values):
    # Initialize the result array
    result1 = np.zeros(3 * 85 * 2 + 2, dtype=np.int32)

    # Handle I values
    for i in range(len(i_values)):
        if i >= 85:  # Ensure we don't exceed the expected 85 values
            break

        # Convert to 48-bit binary representation
        i_val = int(i_values[i])
        if i_val < 0:
            binary_48bit = bin(i_val + 2**48)[2:].zfill(48)
        else:
            binary_48bit = bin(i_val)[2:].zfill(48)

        # Split and store the three 16-bit parts
        result1[i * 3 + 1] = int(binary_48bit[0:16], 2)
        result1[i * 3 + 2] = int(binary_48bit[16:32], 2)
        result1[i * 3 + 3] = int(binary_48bit[32:48], 2)

        # Apply negative value adjustment
        if result1[i * 3 + 1] >= 32768:
            result1[i * 3 + 1] -= 65536
        if result1[i * 3 + 2] >= 32768:
            result1[i * 3 + 2] -= 65536
        if result1[i * 3 + 3] >= 32768:
            result1[i * 3 + 3] -= 65536

    # Handle Q values
    for i in range(len(q_values)):
        if i >= 85:  # Ensure we don't exceed the expected 85 values
            break

        # Convert to 48-bit binary representation
        q_val = int(q_values[i])
        if q_val < 0:
            binary_48bit = bin(q_val + 2**48)[2:].zfill(48)
        else:
            binary_48bit = bin(q_val)[2:].zfill(48)

        # Split and store the three 16-bit parts
        result1[i * 3 + 257] = int(binary_48bit[0:16], 2)
        result1[i * 3 + 258] = int(binary_48bit[16:32], 2)
        result1[i * 3 + 259] = int(binary_48bit[32:48], 2)

        # Apply negative value adjustment
        if result1[i * 3 + 257] >= 32768:
            result1[i * 3 + 257] -= 65536
        if result1[i * 3 + 258] >= 32768:
            result1[i * 3 + 258] -= 65536
        if result1[i * 3 + 259] >= 32768:
            result1[i * 3 + 259] -= 65536

    return result1


def start_generate_raw_data():
    # Read I/Q values from the CSV
    csv_path = "measurement_results_new14.csv"
    i_arrays, q_arrays = read_iq_from_csv(csv_path)

    # Scale the values up by a super large number to make sure they are integers, because fortyeightbit_change_o_single only works with integers
    scale_factor = 1e12

    # Create output directory if it doesn't exist
    os.makedirs("test_output", exist_ok=True)

    # Process each of the 30 I/Q pairs
    for pair_idx in range(len(i_arrays)):
        i_values = i_arrays[pair_idx]
        q_values = q_arrays[pair_idx]

        i_scaled = (i_values * scale_factor).astype(np.int64)
        q_scaled = (q_values * scale_factor).astype(np.int64)

        # Calculate how many groups of 85 rows we have
        total_rows = len(i_scaled)
        num_groups = int(np.ceil(total_rows / 85))

        print(
            f"Pair {pair_idx+1}: Found {total_rows} rows, processing {num_groups} groups of 85 rows each"
        )

        # Process each group
        for group_idx in range(num_groups):
            start_idx = group_idx * 85
            end_idx = min(start_idx + 85, total_rows)

            # Extract current group
            i_group = i_scaled[start_idx:end_idx]
            q_group = q_scaled[start_idx:end_idx]

            # Pad if needed to reach exactly 85 elements
            if len(i_group) < 85:
                i_group = np.pad(i_group, (0, 85 - len(i_group)), "constant")
                q_group = np.pad(q_group, (0, 85 - len(q_group)), "constant")

            # Generate raw data for this group
            raw_data = generate_raw_data_np(i_group, q_group)

            # Verify the conversion for the first few elements (only for the first pair and group)
            if pair_idx == 0 and group_idx == 0:
                ampl, phase, i_verify, q_verify = fortyeightbit_change_o_single(
                    raw_data, 0, 0
                )

                print(
                    "Verification (comparing original to reconstructed, scaled values):"
                )
                for i in range(min(5, len(i_group))):
                    print(f"Original I: {i_group[i]}, Reconstructed I: {i_verify[i]}")
                    print(f"Original Q: {q_group[i]}, Reconstructed Q: {q_verify[i]}")

                    if i_group[i] != i_verify[i] or q_group[i] != q_verify[i]:
                        print(
                            f"{RED}=========== Verification Failed!!!! ================{RESET}"
                        )
                        return

                print(f"{GREEN}=========== Verification Passed ================{RESET}")

            # Determine which qubit this data belongs to and generate appropriate filename
            if pair_idx < 15:
                # First 15 pairs belong to q1
                raw_data_path = f"test_output/raw_data_q1_{group_idx+1}.csv"
            else:
                # Next 15 pairs belong to q2
                raw_data_path = f"test_output/raw_data_q2_{group_idx+1}.csv"

            # Create a DataFrame with a single row where each column is a value
            raw_data_dict = {f"value_{i}": raw_data[i] for i in range(len(raw_data))}
            raw_df = pd.DataFrame([raw_data_dict])
            raw_df.to_csv(raw_data_path, index=False)
            print(f"Raw data saved to {raw_data_path}")


if __name__ == "__main__":
    start_generate_raw_data()
