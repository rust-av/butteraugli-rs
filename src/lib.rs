mod butteraugli;

use crate::butteraugli::open_image;
use std::path::PathBuf;
use structopt::StructOpt;

#[derive(Debug, StructOpt)]
#[structopt(
    name = "butteraugli-rs",
    about = "Rust port of the Google butteraugli metric"
)]
struct CommandLineArgs {
    /// Usage: {*.(png|jpg|jpeg)} {*.(png|jpg|jpeg)}
    #[structopt(parse(from_os_str), required = true, number_of_values = 2)]
    input: Vec<PathBuf>,

    /// In PNM image format {*.(pnm)}
    #[structopt(parse(from_os_str), short, long)]
    output: Option<PathBuf>,
}

pub fn run() {
    dbg!(CommandLineArgs::from_args());

    let k = open_image(CommandLineArgs::from_args().input[1].to_owned());
}
