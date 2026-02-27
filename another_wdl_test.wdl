version 1.0

task run_r {

  input {
    File r_script
  }

  command <<<
    set -euo pipefail
    Rscript ~{r_script}
  >>>

  runtime {
    docker: "rocker/r-base:4.3.1"
  }
}

workflow another_test {

  input {
    File r_script
  }

  call run_r {
    input:
      r_script = r_script
  }
}
