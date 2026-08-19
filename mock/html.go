package mock

import (
	_ "embed"
	"fmt"
	"html/template"
	"os"
	"path/filepath"
	"time"
)

//go:embed report.html.tmpl
var reportTemplate string

// WriteHTML renders a standalone validation report to path.
func WriteHTML(path string, report SuiteReport) error {
	functions := template.FuncMap{
		"duration": func(value time.Duration) string {
			if value < time.Millisecond {
				return fmt.Sprintf("%.0f us", float64(value)/float64(time.Microsecond))
			}
			return fmt.Sprintf("%.2f ms", float64(value)/float64(time.Millisecond))
		},
		"timefmt": func(value time.Time) string { return value.Format("2006-01-02 15:04:05 MST") },
		"status": func(passed bool) string {
			if passed {
				return "通过"
			}
			return "失败"
		},
		"passPercent": func(passed, total int) string {
			if total == 0 {
				return "0"
			}
			return fmt.Sprintf("%.0f", float64(passed)*100/float64(total))
		},
	}
	compiled, err := template.New("validation-report").Funcs(functions).Parse(reportTemplate)
	if err != nil {
		return fmt.Errorf("mock: parse HTML report template: %w", err)
	}
	directory := filepath.Dir(path)
	if err := os.MkdirAll(directory, 0o755); err != nil {
		return fmt.Errorf("mock: create report directory: %w", err)
	}
	file, err := os.Create(path)
	if err != nil {
		return fmt.Errorf("mock: create HTML report: %w", err)
	}
	if err := compiled.Execute(file, report); err != nil {
		_ = file.Close()
		return fmt.Errorf("mock: render HTML report: %w", err)
	}
	if err := file.Close(); err != nil {
		return fmt.Errorf("mock: close HTML report: %w", err)
	}
	return nil
}
