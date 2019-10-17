import pandas as pd
import yaml, sys, os

def create_config(input_csv_fhs, target_size = 2100, index_type = "nextera"):
    configs = []
    for input_csv in input_csv_fhs:
        #read in primer sheet
        df = pd.read_csv(input_csv)
        #filter to defined N7/S5
        df = df[(pd.notna(df["sUMI_Primer_Sequence"])) & (df["Project"] == "IRF3")]
        sub_df = df[["N7_Index","S5_Index","Forward_Primer_2ndRd_Sequence","Reverse_Primer_2ndRd_Sequence","sUMI_Primer_Sequence"]]
        sub_df = sub_df.assign(Run=os.path.basename(input_csv).split(".")[0])
        #iterate
        index_dicts = [sub_df.iloc[i].to_dict() for i in range(sub_df.shape[0])]
        filenames = ["{0}_{1}_{2}".format(df["Sample"].iloc[i],df["N7_Index"].iloc[i],df["S5_Index"].iloc[i]) for i in range(df.shape[0])]
        template_dict = dict(zip(filenames, index_dicts))
        configs.append(template_dict)

    sys.stdout.write(yaml.dump({"datasets": {k: v for d in configs for k, v in d.items()}}))

create_config(sys.argv[1:])
