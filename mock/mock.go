package mock

import (
	"context"
	"flag"
	"fmt"
	"os"
	"time"
)

func Mock() {
	output := flag.String("out", "mock/adjustment-validation.html", "HTML report output path")
	timeout := flag.Duration("timeout", 2*time.Minute, "validation suite timeout")
	flag.Parse()

	ctx, cancel := context.WithTimeout(context.Background(), *timeout)
	defer cancel()
	report := RunAll(ctx)
	if err := WriteHTML(*output, report); err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(2)
	}
	fmt.Printf("validation report: %s\n", *output)
	fmt.Printf("scenarios: %d passed, %d failed, duration %s\n", report.PassCount, report.FailCount, report.Duration.Round(time.Microsecond))
	if !report.Passed {
		for _, scenario := range report.Scenarios {
			if !scenario.Passed {
				fmt.Fprintf(os.Stderr, "failed: %s: %s\n", scenario.ID, scenario.Error)
			}
		}
		os.Exit(1)
	}
}
