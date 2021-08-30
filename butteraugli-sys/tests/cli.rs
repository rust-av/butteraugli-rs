use assert_cmd::prelude::*;
use predicates::prelude::*;
use std::process::Command;

#[test]
fn Run() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("butteraugli-sys")?;

    cmd.arg("test-images/1a.png").arg("test-images/1b.png");
    cmd.assert()
        .success()
        .stdout(predicate::str::contains("133"));

    Ok(())
}
