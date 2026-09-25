import numpy as np


def get_cfus(df, day, condition, culture, strain=None, partner=None, replicate=None):
    subset = df[
        (df["day"] == day) &
        (df["condition"] == condition) &
        (df["culture"] == culture)
    ]

    if strain is not None:
        subset = subset[subset["strain"] == strain]
    if partner is not None:
        subset = subset[subset["partner"] == partner]
    if replicate is not None:
        subset = subset[subset["replicate"] == replicate]

    return subset["CFUs"].values[0]

def get_p_yield(df, condition, partner, replicate):
    p_d0_mono = df[
        (df["day"] == 0) &
        (df["culture"] == "mono") &
        (df["strain"] == partner)
    ]["CFUs"].values[0]

    p_d6_mono = df[
        (df["condition"] == condition) &
        (df["day"] == 6) &
        (df["culture"] == "mono") &
        (df["replicate"] == replicate) &
        (df["strain"] == partner)
    ]["CFUs"].values[0]

    return np.log10(p_d6_mono / p_d0_mono)

def get_oa_yield(df, condition, replicate):
    oa_d0_mono = df[
        (df["day"] == 0) &
        (df["culture"] == "mono") &
        (df["strain"] == "Oa")
    ]["CFUs"].values[0]

    oa_d6_mono = df[
        (df["condition"] == condition) &
        (df["day"] == 6) &
        (df["culture"] == "mono") &
        (df["replicate"] == replicate) &
        (df["strain"] == "Oa")
    ]["CFUs"].values[0]

    return np.log10(oa_d6_mono / oa_d0_mono)

def get_p_on_oa(df, condition, partner, replicate):

    oa_d6_mono = df[
        (df["condition"] == condition) &
        (df["day"] == 6) &
        (df["culture"] == "mono") &
        (df["strain"] == "Oa")
    ]["CFUs"].values[0]

    oa_d6_co = df[
        (df["condition"] == condition) &
        (df["culture"] == "co") &
        (df["replicate"] == replicate) &
        (df["partner"] == partner)
    ]["CFUs"].values[0]

    return np.log10(oa_d6_co / oa_d6_mono)        

def get_oa_on_p(df,condition, partner, replicate):

    p_d6_mono = df[
        (df["condition"] == condition) &
        (df["day"] == 6) &
        (df["culture"] == "mono") &
        (df["strain"] == partner)
    ]["CFUs"].values[0]

    p_d6_co = df[
        (df["condition"] == condition) &
        (df["culture"] == "co") &
        (df["replicate"] == replicate) &
        (df["strain"] == partner)
    ]["CFUs"].values[0]

    return np.log10(p_d6_co / p_d6_mono)

def get_oa_freq(df, condition, partner, replicate):

    oa_d6_co = df[
        (df["condition"] == condition) &
        (df["culture"] == "co") &
        (df["replicate"] == replicate) &
        (df["partner"] == partner)
    ]["CFUs"].values[0]

    p_d6_co = df[
        (df["condition"] == condition) &
        (df["culture"] == "co") &
        (df["replicate"] == replicate) &
        (df["strain"] == partner)
    ]["CFUs"].values[0]

    return oa_d6_co / (oa_d6_co + p_d6_co)

def get_total_cfu(df, condition, partner, replicate):
    
    oa_d6_co = df[
        (df["condition"] == condition) &
        (df["culture"] == "co") &
        (df["replicate"] == replicate) &
        (df["partner"] == partner)
    ]["CFUs"].values[0]

    p_d6_co = df[
        (df["condition"] == condition) &
        (df["culture"] == "co") &
        (df["replicate"] == replicate) &
        (df["strain"] == partner)
    ]["CFUs"].values[0]

    return oa_d6_co + p_d6_co