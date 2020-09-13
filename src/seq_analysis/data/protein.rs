use phf::phf_map;

pub static AMINO_ACID_MASS: phf::Map<u8, f64> = phf_map! {
    b'A' => 89.0932,
    b'C' => 121.1582,
    b'D' => 133.1027,
    b'E' => 147.1293,
    b'F' => 165.1891,
    b'G' => 75.0666,
    b'H' => 155.1546,
    b'I' => 131.1729,
    b'K' => 146.1876,
    b'L' => 131.1729,
    b'M' => 149.2113,
    b'N' => 132.1179,
    b'O' => 255.3134,
    b'P' => 115.1305,
    b'Q' => 146.1445,
    b'R' => 174.201,
    b'S' => 105.0926,
    b'T' => 119.1192,
    b'U' => 168.0532,
    b'V' => 117.1463,
    b'W' => 204.2252,
    b'Y' => 181.1885,
};

pub static AMINO_ACID_MASS_MONOISOTOPIC: phf::Map<u8, f64> = phf_map! {
    b'A'=> 89.047678,
    b'C'=> 121.019749,
    b'D'=> 133.037508,
    b'E'=> 147.053158,
    b'F'=> 165.078979,
    b'G'=> 75.032028,
    b'H'=> 155.069477,
    b'I'=> 131.094629,
    b'K'=> 146.105528,
    b'L'=> 131.094629,
    b'M'=> 149.051049,
    b'N'=> 132.053492,
    b'O'=> 255.158292,
    b'P'=> 115.063329,
    b'Q'=> 146.069142,
    b'R'=> 174.111676,
    b'S'=> 105.042593,
    b'T'=> 119.058243,
    b'U'=> 168.964203,
    b'V'=> 117.078979,
    b'W'=> 204.089878,
    b'Y'=> 181.073893,
};

/// Flexibility
///
/// Normalized flexibility parameters (B-values), average
///
/// Vihinen M., Torkkila E., Riikonen P. Proteins. 19(2):141-9(1994).
pub static AMINO_ACID_FLEX: phf::Map<u8, f32> = phf_map! {
    b'A' => 0.984,
    b'C' => 0.906,
    b'E' => 1.094,
    b'D' => 1.068,
    b'G' => 1.031,
    b'F' => 0.915,
    b'I' => 0.927,
    b'H' => 0.950,
    b'K' => 1.102,
    b'M' => 0.952,
    b'L' => 0.935,
    b'N' => 1.048,
    b'Q' => 1.037,
    b'P' => 1.049,
    b'S' => 1.046,
    b'R' => 1.008,
    b'T' => 0.997,
    b'W' => 0.904,
    b'V' => 0.931,
    b'Y' => 0.929
};

pub const N_TERM_PKA_DEFAULT: f32 = 7.5;
pub const C_TERM_PKA_DEFAULT: f32 = 3.55;

pub enum Charge {
    Positive,
    Negative,
}

pub static PKA: phf::Map<u8, (f32, Charge)> = phf_map! {
    b'K'=> (10.0, Charge::Positive),
    b'R'=> (12.0, Charge::Positive),
    b'H'=> (5.98, Charge::Positive),
    b'D'=> (4.05, Charge::Negative),
    b'E'=> (4.45, Charge::Negative),
    b'C'=> (9.00, Charge::Negative),
    b'Y'=> (10.0, Charge::Negative),
};

pub static PKA_N_TERM: phf::Map<u8, f32> = phf_map! {
    b'A' => 7.59,
    b'M' => 7.00,
    b'S' => 6.93,
    b'P' => 8.36,
    b'T' => 6.82,
    b'V' => 7.44,
    b'E' => 7.70,
};

pub static PKA_C_TERM: phf::Map<u8, f32> = phf_map! {
    b'D' => 4.55,
    b'E' => 4.75,
};
