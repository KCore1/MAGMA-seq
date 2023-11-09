import pandas as pd

input_file = "/projects/brpe7306/fab-library-barcoding/MLE/input/rep1/HA_rep1_collapsed.csv"
output_file = "/projects/brpe7306/fab-library-barcoding/MLE/input/rep1/HA_rep1_collapsed_bright_adjust.csv"

alpha_adjusted = 0.5345
adjusted_labels = [25, 100, 250, 750, 1000]

df = pd.read_csv(input_file)

df['brightness_adjust'] = df['Concentration'].apply(lambda x: 1 if x not in adjusted_labels else alpha_adjusted)
df.to_csv(output_file, index=False)