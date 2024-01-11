

def count_commas_in_csv(file_path):
    with open(file_path, 'r') as file:
        for line_number, line in enumerate(file, start=1):
            comma_count = line.count(',')
            print(f"Line {line_number}: {comma_count} commas")

# Usage example:
file = '8_2_R3e-05_v4e-01_cor0.63_mu0.1_rho2.25_k3e+00_Ha5e-12_dt4e-10_simData.csv'
folder = '/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/error2Test2/N_10/T_10/'
count_commas_in_csv(folder+file)