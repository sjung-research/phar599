// Build the report .docx from its Markdown source. Handles the subset the report uses:
// title (#), headings (##), paragraphs with **bold** / *italic*, bullet and numbered
// lists, pipe tables, ![alt](path.png) images, and italic caption lines.
// Usage: node reports/build_report_docx.js reports/<report>.md reports/<report>.docx
const fs = require("fs");
const path = require("path");
const {
  Document, Packer, Paragraph, TextRun, HeadingLevel, Table, TableRow, TableCell,
  WidthType, ShadingType, BorderStyle, ImageRun, AlignmentType, LevelFormat,
  Footer, PageNumber, LineRuleType,
} = require("docx");

const [, , srcArg, outArg] = process.argv;
const src = path.resolve(srcArg);
const out = path.resolve(outArg);
const baseDir = path.dirname(src);
const lines = fs.readFileSync(src, "utf8").replace(/\r/g, "").split("\n");

const FONT = "Arial";
const TEXT_W = 9360; // 6.5 in in DXA (US Letter, 1 in margins)
const IMG_W_PX = 624; // 6.5 in at 96 dpi
const BORDER = { style: BorderStyle.SINGLE, size: 4, color: "BFBFBF" };

function runs(text, base = {}) {
  // **bold** and *italic*; underscores stay literal (e.g. _Lng)
  const parts = [];
  const re = /(\*\*[^*]+\*\*|\*[^*]+\*)/g;
  let last = 0, m;
  while ((m = re.exec(text)) !== null) {
    if (m.index > last) parts.push(new TextRun({ text: text.slice(last, m.index), ...base }));
    const t = m[0];
    parts.push(t.startsWith("**")
      ? new TextRun({ text: t.slice(2, -2), ...base, bold: true })
      : new TextRun({ text: t.slice(1, -1), ...base, italics: true }));
    last = m.index + t.length;
  }
  if (last < text.length) parts.push(new TextRun({ text: text.slice(last), ...base }));
  return parts;
}

function pngSize(file) {
  const b = fs.readFileSync(file);
  if (b.toString("ascii", 1, 4) !== "PNG") throw new Error(`not a PNG: ${file}`);
  return { w: b.readUInt32BE(16), h: b.readUInt32BE(20), data: b };
}

function table(rows) {
  const cells = rows.map((r) => r.replace(/^\||\|$/g, "").split("|").map((c) => c.trim()));
  const header = cells[0];
  const body = cells.slice(2); // skip the --- row
  const n = header.length;
  const col = (i) => [header, ...body].map((r) => r[i] || "");
  const longestWord = (i) => Math.max(...col(i).flatMap((c) => c.split(/\s+/)).map((w) => w.length));
  const weights = header.map((_, i) =>
    Math.min(40, Math.max(longestWord(i) + 2, ...col(i).map((c) => c.length))));
  const total = weights.reduce((a, b) => a + b, 0);
  const widths = weights.map((w) => Math.floor((w / total) * TEXT_W));
  widths[n - 1] += TEXT_W - widths.reduce((a, b) => a + b, 0);
  const mk = (r, isHeader) => new TableRow({
    tableHeader: isHeader,
    cantSplit: true,
    children: r.map((c, i) => new TableCell({
      width: { size: widths[i], type: WidthType.DXA },
      borders: { top: BORDER, bottom: BORDER, left: BORDER, right: BORDER },
      shading: isHeader ? { type: ShadingType.CLEAR, color: "auto", fill: "E8EEF5" } : undefined,
      margins: { top: 60, bottom: 60, left: 90, right: 90 },
      children: [new Paragraph({ spacing: { after: 0 }, children: runs(c, isHeader ? { size: 17, bold: true } : { size: 17 }) })],
    })),
  });
  return new Table({
    width: { size: TEXT_W, type: WidthType.DXA },
    columnWidths: widths,
    rows: [mk(header, true), ...body.map((r) => mk(r, false))],
  });
}

const numbering = [{
  reference: "bullets",
  levels: [{ level: 0, format: LevelFormat.BULLET, text: "•", alignment: AlignmentType.LEFT,
             style: { paragraph: { indent: { left: 540, hanging: 270 } } } }],
}];
let listNo = 0;
const children = [];
let i = 0;
while (i < lines.length) {
  const line = lines[i];
  if (!line.trim()) { i++; continue; }
  if (line.startsWith("# ")) {
    children.push(new Paragraph({ heading: HeadingLevel.TITLE, children: runs(line.slice(2)) }));
  } else if (line.startsWith("## ")) {
    children.push(new Paragraph({ heading: HeadingLevel.HEADING_1, children: runs(line.slice(3)) }));
  } else if (line.startsWith("|")) {
    const rows = [];
    while (i < lines.length && lines[i].startsWith("|")) rows.push(lines[i++]);
    children.push(table(rows));
    children.push(new Paragraph({ spacing: { after: 60 }, children: [] }));
    continue;
  } else if (/^!\[/.test(line)) {
    const m = line.match(/^!\[([^\]]*)\]\(([^)]+)\)/);
    const img = pngSize(path.join(baseDir, m[2]));
    const h = Math.round((img.h / img.w) * IMG_W_PX);
    children.push(new Paragraph({
      alignment: AlignmentType.CENTER, keepNext: true,
      spacing: { before: 120, after: 60, line: 240, lineRule: LineRuleType.AUTO },
      children: [new ImageRun({ type: "png", data: img.data, transformation: { width: IMG_W_PX, height: h },
                                altText: { title: m[1], description: m[1], name: path.basename(m[2]) } })],
    }));
  } else if (line.startsWith("- ")) {
    while (i < lines.length && lines[i].startsWith("- ")) {
      children.push(new Paragraph({ numbering: { reference: "bullets", level: 0 }, children: runs(lines[i].slice(2)) }));
      i++;
    }
    continue;
  } else if (/^\d+\. /.test(line)) {
    const ref = `numbers-${listNo++}`; // each list restarts at 1
    numbering.push({
      reference: ref,
      levels: [{ level: 0, format: LevelFormat.DECIMAL, text: "%1.", alignment: AlignmentType.LEFT,
                 style: { paragraph: { indent: { left: 540, hanging: 300 } } } }],
    });
    while (i < lines.length && /^\d+\. /.test(lines[i])) {
      children.push(new Paragraph({ numbering: { reference: ref, level: 0 }, children: runs(lines[i].replace(/^\d+\. /, "")) }));
      i++;
    }
    continue;
  } else if (/^\*[^*].*\*$/.test(line.trim())) {
    // figure caption
    children.push(new Paragraph({ spacing: { after: 220 }, children: runs(line.trim(), { size: 18, color: "444444" }) }));
  } else {
    children.push(new Paragraph({ children: runs(line) }));
  }
  i++;
}

const doc = new Document({
  creator: "Sewon",
  title: "Mouse anti-GFAP concentration curve: first quantitative report",
  styles: {
    default: { document: { run: { font: FONT, size: 21 }, paragraph: { spacing: { after: 120, line: 276, lineRule: LineRuleType.AUTO } } } },
    paragraphStyles: [
      { id: "Title", name: "Title", basedOn: "Normal", next: "Normal",
        run: { font: FONT, size: 34, bold: true, color: "0D366B" },
        paragraph: { spacing: { before: 0, after: 160 }, keepNext: true } },
      { id: "Heading1", name: "Heading 1", basedOn: "Normal", next: "Normal", quickFormat: true,
        run: { font: FONT, size: 26, bold: true, color: "0D366B" },
        paragraph: { spacing: { before: 280, after: 120 }, keepNext: true, outlineLevel: 0 } },
    ],
  },
  numbering: { config: numbering },
  sections: [{
    properties: { page: { size: { width: 12240, height: 15840 },
                          margin: { top: 1440, bottom: 1440, left: 1440, right: 1440 } } },
    footers: { default: new Footer({ children: [new Paragraph({ alignment: AlignmentType.CENTER,
      children: [new TextRun({ children: ["Page ", PageNumber.CURRENT, " of ", PageNumber.TOTAL_PAGES], size: 16, color: "777777" })] })] }) },
    children,
  }],
});

Packer.toBuffer(doc).then((buf) => { fs.writeFileSync(out, buf); console.log(`wrote ${out} (${buf.length} bytes)`); });
