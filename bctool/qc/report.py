"""
QC Report generator using HTML+CSS bars (no Unicode/encoding issues).
"""

from pathlib import Path


class QCReport:
    def __init__(self, analyzer, output_dir):
        self.metrics = analyzer.metrics
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.analyzer = analyzer

    def generate_html(self, filename="qc_report.html"):
        m = self.metrics
        grade = self.analyzer.get_grade()
        gclass = "pass" if "PASS" in grade else ("warn" if "WARNING" in grade else "fail")

        parts = []
        parts.append("""<!DOCTYPE html>
<html><head><meta charset="utf-8">
<title>QC Report</title>
<style>
body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;
  max-width:1100px;margin:0 auto;padding:20px;background:#f5f5f5;color:#333;}
h1{color:#333;border-bottom:2px solid #4a90d9;padding-bottom:10px;}
.module{background:#fff;border-radius:8px;box-shadow:0 2px 4px rgba(0,0,0,0.1);
  margin:15px 0;padding:20px;overflow-x:auto;}
table{border-collapse:collapse;width:100%;}
th{background:#4a90d9;color:#fff;padding:8px 12px;text-align:left;}
td{padding:6px 12px;border-bottom:1px solid #eee;}
.bar-wrap{display:flex;align-items:center;margin:2px 0;}
.bar-label{width:60px;font-size:12px;color:#666;text-align:right;padding-right:8px;}
.bar-fill{height:16px;border-radius:3px;min-width:2px;}
.bar-a{background:#4CAF50;} .bar-c{background:#2196F3;}
.bar-g{background:#FF9800;} .bar-t{background:#f44336;}
.bar-n{background:#9E9E9E;} .bar-q{background:#4a90d9;}
.bar-gc{background:#43a047;}
.grade{display:inline-block;padding:8px 20px;border-radius:20px;
  font-weight:bold;font-size:18px;}
.pass{background:#e8f5e9;color:#2e7d32;}
.warn{background:#fff3e0;color:#e65100;}
.fail{background:#ffebee;color:#c62828;}
.warn-box{background:#fff3e0;border-left:4px solid #ff9800;padding:10px 15px;margin:10px 0;border-radius:4px;}
.fail-box{background:#ffebee;border-left:4px solid #f44336;padding:10px 15px;margin:10px 0;border-radius:4px;}
</style></head><body>
<h1>QC Report</h1>
<div class="module">
  <span class="grade """ + gclass + """">""" + grade + """</span>
  <span style="margin-left:15px;">Reads: """ + str(m['total_reads']) + """ | Bases: """ + f"{m['total_bases']:,}" + """</span>
</div>""")

        parts.append(self._section_summary())
        parts.append(self._section_per_base_qual())
        parts.append(self._section_per_base_content())
        parts.append(self._section_gc_content())
        parts.append(self._section_per_seq_qual())
        parts.append(self._section_length_dist())
        parts.append(self._section_overrepresented())
        parts.append(self._section_adapters())
        parts.append("</body></html>")

        html = "\n".join(parts)
        path = self.output_dir / filename
        with open(path, "w", encoding="utf-8") as f:
            f.write(html)
        print(f"  [QC] Report: {path}")
        return str(path)

    def _section_summary(self):
        m = self.metrics
        bp = m.get("base_pct", {})
        qs = m.get("qual_summary", {})
        ls = m.get("length_summary", {})
        gc = m.get("gc_summary", {})
        h = '<div class="module"><h2>Summary</h2><table>'
        h += f"<tr><td>Total reads</td><td>{m['total_reads']}</td></tr>"
        h += f"<tr><td>Total bases</td><td>{m['total_bases']:,}</td></tr>"
        h += f"<tr><td>Read length</td><td>{ls.get('mean',0):.1f} (N50:{ls.get('n50',0)})</td></tr>"
        h += f"<tr><td>GC content</td><td>{gc.get('mean',0):.1f}%</td></tr>"
        h += f"<tr><td>Mean quality</td><td>{qs.get('mean',0):.1f}</td></tr>"
        for b in "ACGTN":
            h += f"<tr><td>{b}</td><td>{bp.get(b,0):.2f}%</td></tr>"
        for ad in m.get("detected_adapters", []):
            h += f'<tr class="warn-box"><td>Adapter</td><td>{ad["adapter"]} ({ad["pct"]:.1f}%)</td></tr>'
        h += "</table></div>"
        return h

    def _section_per_base_qual(self):
        sq = self.metrics.get("per_base_qual_summary", {})
        if not sq:
            return ""
        max_pos = max(sq.keys()) + 1
        h = '<div class="module"><h2>Per-Base Quality</h2>'
        h += '<div class="bar-wrap"><span class="bar-label"></span>'
        for i in range(max_pos):
            h += f'<span style="width:12px;font-size:10px;text-align:center;">{i%10}</span>'
        h += '</div>'
        for q in [40, 35, 30, 25, 20, 15, 10, 5, 0]:
            label = f"{q:>2}-{q+4}" if q > 0 else " 0-4"
            h += f'<div class="bar-wrap"><span class="bar-label">{label}</span>'
            for pos in range(max_pos):
                if pos in sq:
                    s = sq[pos]
                    w = 12
                    color = "#43a047" if s["mean"] >= 28 else ("#ff9800" if s["mean"] >= 20 else "#f44336")
                    h += f'<span style="width:{w}px;height:14px;display:inline-block;background:{color};border-radius:1px;margin:0 0;"></span>'
                else:
                    h += f'<span style="width:12px;display:inline-block;"></span>'
            h += '</div>'
        means = [s["mean"] for p,s in sq.items() if "mean" in s]
        maxes = [s["max"] for p,s in sq.items() if "max" in s]
        if means:
            h += f'<p>Mean range: {min(means):.0f}-{max(means):.0f} | Max: {max(maxes):.0f}</p>'
        h += '</div>'
        return h

    def _section_per_base_content(self):
        content = self.metrics.get("per_base_content", {})
        if not content:
            return ""
        max_pos = max(content.keys()) + 1
        h = '<div class="module"><h2>Per-Base Sequence Content</h2>'
        for base, color, label in [("A","#4CAF50","A"),("T","#f44336","T")]:
            h += f'<div class="bar-wrap"><span class="bar-label">{label}</span>'
            for pos in range(max_pos):
                if pos in content:
                    c = content[pos]
                    total = sum(c.values()) or 1
                    pct = c[base] / total * 100
                    h += f'<span style="width:12px;height:{pct*0.4:.0f}px;display:inline-block;background:{color};border-radius:1px;vertical-align:bottom;"></span>'
                else:
                    h += '<span style="width:12px;display:inline-block;"></span>'
            h += '</div>'
        h += '<p style="font-size:12px;color:#666;">Green=A, Red=T (G/C similar)</p></div>'
        return h

    def _section_gc_content(self):
        gc = self.metrics.get("gc_bins", {})
        if not gc:
            return ""
        maxc = max(gc.values()) or 1
        scale = 200 / maxc
        gc_mean = self.metrics.get("gc_summary", {}).get("mean", 50)
        h = '<div class="module"><h2>GC Distribution</h2>'
        for g in range(0, 101, 5):
            count = sum(gc.get(i,0) for i in range(g, g+5))
            w = int(count * scale)
            h += f'<div class="bar-wrap"><span class="bar-label">{g:>2}%</span>'
            h += f'<span class="bar-fill bar-gc" style="width:{min(w,400)}px;"></span>'
            h += f'<span style="margin-left:5px;font-size:11px;color:#666;">{count}</span></div>'
        h += f'<p>Mean GC: {gc_mean:.1f}%</p></div>'
        return h

    def _section_per_seq_qual(self):
        qd = self.metrics.get("qual_bins", {})
        if not qd:
            return ""
        maxc = max(qd.values()) or 1
        scale = 200 / maxc
        h = '<div class="module"><h2>Per-Sequence Quality</h2>'
        for q in range(0, 42, 2):
            count = sum(qd.get(i,0) for i in range(max(0,q-1), min(42,q+2)))
            w = int(count * scale)
            h += f'<div class="bar-wrap"><span class="bar-label">{q:>2}</span>'
            h += f'<span class="bar-fill bar-q" style="width:{min(w,400)}px;"></span>'
            h += f'<span style="margin-left:5px;font-size:11px;color:#666;">{count}</span></div>'
        h += '</div>'
        return h

    def _section_length_dist(self):
        lengths = self.metrics.get("read_lengths", {})
        if not lengths:
            return ""
        maxc = max(lengths.values()) or 1
        scale = 200 / maxc
        ls = self.metrics.get("length_summary", {})
        h = '<div class="module"><h2>Read Length Distribution</h2>'
        for length in sorted(lengths):
            count = lengths[length]
            w = int(count * scale)
            h += f'<div class="bar-wrap"><span class="bar-label">{length}</span>'
            h += f'<span class="bar-fill bar-q" style="width:{min(w,400)}px;"></span>'
            h += f'<span style="margin-left:5px;font-size:11px;color:#666;">{count}</span></div>'
        h += f'<p>Min:{ls.get("min",0)} Max:{ls.get("max",0)} Mean:{ls.get("mean",0):.1f} N50:{ls.get("n50",0)}</p></div>'
        return h

    def _section_overrepresented(self):
        over = self.metrics.get("overrepresented", [])
        if not over:
            return ""
        h = '<div class="module"><h2>Overrepresented Sequences</h2><div class="warn-box">Found overrepresented sequences.</div><table><tr><th>Count</th><th>%</th><th>Sequence</th></tr>'
        for item in over[:20]:
            h += f"<tr><td>{item['count']}</td><td>{item['pct']:.2f}%</td><td style='font-family:monospace;font-size:11px;'>{item['sequence']}</td></tr>"
        h += "</table></div>"
        return h

    def _section_adapters(self):
        ad = self.metrics.get("detected_adapters", [])
        if not ad:
            return ""
        h = '<div class="module"><h2>Adapter Contamination</h2><div class="fail-box">Adapter detected! Use --trim-adapters.</div><table><tr><th>Adapter</th><th>%</th></tr>'
        for a in ad:
            h += f"<tr><td>{a['adapter']}</td><td>{a['pct']:.2f}%</td></tr>"
        h += "</table></div>"
        return h
