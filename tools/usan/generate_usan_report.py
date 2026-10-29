#!/usr/bin/env python3
"""Generate the USAN / WHO INN antibody evaluation standard guide (Word).

The guide is the project's written standard for (a) how humanization degree is
estimated from sequence data, (b) how USAN / WHO INN antibody nomenclature is
interpreted, and (c) how to phrase results in papers and reports. It is a
documentation artifact -- it does not run the pipeline.

Facts are aligned with:
  * WHO INN Working Document 22-542 (May 2022), the new mAb nomenclature
    scheme that retired the single "-mab" stem in favour of four format-based
    stems (-tug / -bart / -ment / -mig);
  * the legacy origin-based suffixes -umab / -zumab / -ximab / -omab;
  * IMGT unique numbering (Lefranc 2003), as implemented in
    scripts/humanize/imgt_numbering.py;
  * the academic (non-official) FR-identity conventions 95 / 85 / 70 % used by
    tools/humanness_score/evaluate_humanness.py.

Usage:
    python3 tools/usan/generate_usan_report.py [--output PATH] [--date YYYY-MM-DD]
"""

from __future__ import annotations

import argparse
import datetime
import os
import sys

try:
    from docx import Document
    from docx.enum.table import WD_TABLE_ALIGNMENT
    from docx.enum.text import WD_ALIGN_PARAGRAPH
    from docx.oxml import OxmlElement
    from docx.oxml.ns import qn
    from docx.shared import Cm, Pt, RGBColor
except ImportError:  # pragma: no cover - optional dependency
    sys.stderr.write(
        "ERROR: python-docx is required to generate the USAN guide.\n"
        "       install it with:  pip install python-docx\n"
    )
    raise SystemExit(2)


DOC_VERSION = "v2.0"
DOC_ID = "USAN-GUIDE-001"

ASCII_FONT_BODY = "Times New Roman"
ASCII_FONT_HEAD = "Arial"
CJK_FONT_BODY = "宋体"
CJK_FONT_HEAD = "黑体"

HEADER_FILL = "D9E2F3"
META_FILL = "F2F2F2"

# ---------------------------------------------------------------------------
# Table data (single source of truth for every table in the guide)
# ---------------------------------------------------------------------------

GLOSSARY = [
    ("USAN", "United States Adopted Name；美国采用名称，由 AMA 药品命名委员会（历史上与 USP 共同）审定。"),
    ("INN", "International Nonproprietary Names；WHO 主持的国际非专有名体系。"),
    ("AMA / USP", "American Medical Association / United States Pharmacopeia（USAN 历史共同管理方）。"),
    ("IMGT", "IMMunoGeneTics 数据库及其唯一编号体系（Lefranc 2003）。"),
    ("Kabat", "Kabat 等人提出的抗体编号体系（1979 / 1991）。"),
    ("FR / CDR", "框架区（framework region）/ 互补决定区（complementarity-determining region）。"),
    ("FR4", "V 区 CDR3 之后的框架段，由 J 基因编码（本管线评估时排除）。"),
    ("VH / VL / VHH", "重链可变区 / 轻链可变区 / 重链单域可变区（纳米抗体）。"),
    ("Germline", "种系 V 基因；与查询序列比对计算同源性。"),
    ("CDR grafting", "CDR 移植：将供体 CDR 挂载到人源框架上的人源化方法。"),
    ("Back-mutation", "回复突变：将人源化后丢失的关键供体残基恢复。"),
    ("CVI", "Canonical + Vernier + Interface 残基集合，用于结构保守性同源性评估。"),
]

SUFFIX_TABLE = [
    ("-umab", "全人源（fully human）", "转基因小鼠 / 噬菌体展示等全人源文库",
     "Adalimumab、Nivolumab、Denosumab"),
    ("-zumab", "人源化（humanized）", "CDR 移植（供体 CDR + 人源框架）",
     "Trastuzumab、Bevacizumab、Pembrolizumab"),
    ("-ximab", "嵌合（chimeric）", "鼠源 V 区 + 人源恒定区",
     "Rituximab、Infliximab"),
    ("-omab", "鼠源（murine）", "杂交瘤等传统技术，全鼠抗体",
     "Muromonab-CD3"),
]

NEW_STEMS_TABLE = [
    ("-tug", "unmodified immunoglobulins",
     "单特异性全长免疫球蛋白，恒定区与天然序列一致（允许工程化聚糖、C 端 Lys 缺失）；"
     "涵盖天然、嵌合、人源化、全人源及骆驼科重链抗体。",
     "天然 / 嵌合 / 人源化 / 全人源全长抗体"),
    ("-bart", "artificial immunoglobulins",
     "单特异性全长免疫球蛋白，恒定区（含铰链）存在工程化氨基酸改变，"
     "如 IgG4 S-P 铰链、Fc 效应功能削弱、FcRn 结合改造、新糖基化位点等。",
     "Fc 工程化 / 恒定区改造抗体"),
    ("-ment", "immunoglobulin fragments",
     "不属于 -tug / -bart 的单特异性片段，至少含一个参与结合的可变结构域，"
     "恒定区可完整、部分或缺失（如 scFv、Fab、无 Fc 的 scFv-Fc 构建体）。",
     "抗体片段 / 无 Fc 构建体"),
    ("-mig", "multi-specific immunoglobulins",
     "双特异性及多特异性免疫球蛋白，不论格式（常规或工程化）、类型（全长或片段）"
     "或形态（是否带延伸）；不包含仅靠单一 CDR 组交叉反应的抗体。",
     "双特异性 / 多特异性抗体"),
]

INFIX_TABLE = [
    ("-ta-", "tumour", "肿瘤"), ("-ci-", "cardiovascular", "心血管"),
    ("-de-", "metabolic / endocrine", "代谢 / 内分泌"), ("-eni-", "enzyme inhibition", "酶抑制"),
    ("-gro-", "growth factor (receptor)", "生长因子（受体）"), ("-ki-", "cytokine (receptor)", "细胞因子（受体）"),
    ("-ne-", "neural", "神经系统"), ("-os-", "bone", "骨"),
    ("-vi-", "viral", "病毒"), ("-ba-", "bacterial", "细菌"),
    ("-fung-", "fungal", "真菌"), ("-ler-", "allergen", "过敏原"),
    ("-sto-", "immunostimulatory", "免疫刺激"), ("-pru-", "immunosuppressive", "免疫抑制"),
    ("-toxa-", "toxin", "毒素"), ("-ami-", "SAP / amyloidosis", "血清淀粉样蛋白 / 淀粉样变"),
    ("-vet-", "veterinary use", "兽用"),
]

CDR_TABLE = [
    ("CDR-H1", "27-38", "31-35", "IMGT ✓"),
    ("CDR-H2", "56-65", "50-65", "IMGT ✓"),
    ("CDR-H3", "105-117", "95-102", "IMGT ✓"),
    ("CDR-L1", "27-38", "24-34", "IMGT ✓"),
    ("CDR-L2", "56-65", "50-56", "IMGT ✓"),
    ("CDR-L3", "105-117", "89-97", "IMGT ✓"),
]

IMGT_REGION_TABLE = [
    ("FR1", "1-26", "CDR1", "27-38", "FR2", "39-55"),
    ("CDR2", "56-65", "FR3", "66-104", "CDR3", "105-117"),
    ("FR4", "118-128", "", "", "", ""),
]

THRESHOLD_TABLE = [
    ("FR 同源性 ≥ 95%", "likely fully-human-like", "全人源（-umab）",
     "转基因小鼠 / 全人源文库筛选"),
    ("85% ≤ FR 同源性 < 95%", "likely humanized-like", "人源化（-zumab）",
     "CDR 移植 + 回复突变"),
    ("70% ≤ FR 同源性 < 85%", "likely chimeric-like", "嵌合（-ximab）",
     "鼠 V 区 + 人 C 区"),
    ("FR 同源性 < 70%", "likely murine-like", "鼠源（-omab）",
     "全鼠抗体"),
]

FULLY_HUMAN_EXAMPLES = [
    ("Adalimumab", "TNF-α", "类风湿关节炎等", "2002", "噬菌体展示"),
    ("Nivolumab", "PD-1", "黑色素瘤等", "2014", "转基因小鼠"),
    ("Denosumab", "RANKL", "骨质疏松 / 骨转移", "2010", "转基因小鼠"),
    ("Avelumab", "PD-L1", "Merkel 细胞癌等", "2017", "转基因小鼠"),
]

HUMANIZED_EXAMPLES = [
    ("Trastuzumab", "HER2", "乳腺癌", "1998", "CDR 移植（鼠源 4D5）"),
    ("Bevacizumab", "VEGF-A", "结直肠癌等", "2004", "CDR 移植（鼠源 A4.6.1）"),
    ("Natalizumab", "α4 整合素", "多发性硬化", "2004", "CDR 移植"),
    ("Pembrolizumab", "PD-1", "黑色素瘤等", "2014", "人源化（平台信息见原研方披露）"),
]

JSON_MAPPING_TABLE = [
    ("germline_imgt.fr_identity", "FR 同源性（IMGT，FR1+FR2+FR3）", "§3.4"),
    ("germline_imgt.cdr_identity", "CDR 同源性（IMGT，CDR1+CDR2，不含 CDR3）", "§3.4"),
    ("germline_imgt.all_identity", "合并同源性（FR1+CDR1+FR2+CDR2+FR3）", "§3.4"),
    ("germline_imgt.usan_naming_class", "学术估计标签（JSON 简化为 85 / 70 两档）", "§3.5"),
    ("cvi_homology", "CVI 残基同源性（BI 2024；本管线 config.py 位置集）", "§5.3"),
    ("backmutations / variants", "回复突变位点与 V0-V3 变体序列", "§1.1"),
    ("tools/humanness_score 输出", "同源性 + 四档学术分级（95 / 85 / 70）", "§3.5"),
]

REFERENCES = [
    "World Health Organization. International Nonproprietary Names (INN) for "
    "Biological Products. WHO, Geneva.",
    "World Health Organization. New INN monoclonal antibody (mAb) nomenclature "
    "scheme. INN Working Document 22-542, May 2022.",
    "Lefranc M-P. IMGT unique numbering for immunoglobulin and T cell receptor "
    "variable domains and Ig superfamily V-like domains. Dev Comp Immunol. "
    "2003;27(1):55-77.",
    "Lefranc M-P, et al. IMGT, the international ImMunoGeneTics database. "
    "Nucleic Acids Res. 1999;27(1):209-212.",
    "Kabat EA, Wu TT, Perry HM, Gottesman KS, Foeller C. Sequences of Proteins "
    "of Immunological Interest. 5th ed. NIH Publication No. 91-3242, 1991.",
    "Chothia C, Lesk AM, Tramontano A, et al. Conformations of immunoglobulin "
    "hypervariable regions. Nature. 1989;342(6245):877-883.",
    "Chothia C, Lesk AM. Canonical structures for the hypervariable regions of "
    "immunoglobulins. J Mol Biol. 1987;196(4):901-917.",
    "Al-Lazikani B, Lesk AM, Chothia C. Standard conformations for the canonical "
    "structures of immunoglobulins. J Mol Biol. 1997;273(4):927-948.",
    "Foote J, Winter G. Antibody framework residues affecting the conformation "
    "of the hypervariable loops. J Mol Biol. 1992;224(2):487-499.",
    "Chothia C, et al. Residues in the VH/VL packing interface. J Mol Biol. "
    "1985;186(3):651-663.",
    "Vargas-Madrazo E, Paz-Garcia E. An extended structural analysis of "
    "antibody-antigen complexes. J Mol Biol. 2003;330(4):783-798.",
    "Meaney A, et al. Vernier, canonical and interface residues are functionally "
    "constrained. MAbs. 2024.",
    "Boehringer Ingelheim. Framework-matrix transplantation and the CVI "
    "homology metric. J Biol Chem. 2024;300:105555.",
]


# ---------------------------------------------------------------------------
# Word helpers
# ---------------------------------------------------------------------------

def _style_run(run, *, size=None, bold=None, color=None,
               east_asia=CJK_FONT_BODY, ascii_font=ASCII_FONT_BODY):
    run.font.name = ascii_font
    if size is not None:
        run.font.size = Pt(size)
    if bold is not None:
        run.font.bold = bold
    if color is not None:
        run.font.color.rgb = color
    rpr = run._element.get_or_add_rPr()
    rfonts = rpr.get_or_add_rFonts()
    rfonts.set(qn("w:eastAsia"), east_asia)
    return run


def _set_style_fonts(style, ascii_font, east_asia, size=None):
    style.font.name = ascii_font
    if size is not None:
        style.font.size = Pt(size)
    rpr = style.element.get_or_add_rPr()
    rfonts = rpr.get_or_add_rFonts()
    rfonts.set(qn("w:eastAsia"), east_asia)


def _shade_cell(cell, fill):
    tc_pr = cell._tc.get_or_add_tcPr()
    shd = OxmlElement("w:shd")
    shd.set(qn("w:val"), "clear")
    shd.set(qn("w:color"), "auto")
    shd.set(qn("w:fill"), fill)
    tc_pr.append(shd)


def _repeat_header_row(table):
    tr = table.rows[0]._tr
    tr_pr = tr.get_or_add_trPr()
    header = OxmlElement("w:tblHeader")
    header.set(qn("w:val"), "true")
    tr_pr.append(header)


def _add_field(paragraph, instr):
    """Insert a Word field (e.g. PAGE, TOC) into a paragraph."""
    run = paragraph.add_run()
    begin = OxmlElement("w:fldChar")
    begin.set(qn("w:fldCharType"), "begin")
    run._r.append(begin)
    run = paragraph.add_run()
    text = OxmlElement("w:instrText")
    text.set(qn("xml:space"), "preserve")
    text.text = instr
    run._r.append(text)
    run = paragraph.add_run()
    separate = OxmlElement("w:fldChar")
    separate.set(qn("w:fldCharType"), "separate")
    run._r.append(separate)
    run = paragraph.add_run("更新域后显示")
    _style_run(run, size=10.5, color=RGBColor(0x80, 0x80, 0x80))
    run = paragraph.add_run()
    end = OxmlElement("w:fldChar")
    end.set(qn("w:fldCharType"), "end")
    run._r.append(end)


class GuideBuilder:
    """Small facade over python-docx with captions and consistent styling."""

    def __init__(self, doc):
        self.doc = doc
        self.table_no = 0

    # -- structural elements ------------------------------------------------

    def heading(self, text, level):
        p = self.doc.add_heading("", level=level)
        run = p.add_run(text)
        size = {1: 16, 2: 14, 3: 12}.get(level, 12)
        _style_run(run, size=size, bold=True, east_asia=CJK_FONT_HEAD,
                   ascii_font=ASCII_FONT_HEAD,
                   color=RGBColor(0, 0, 0))
        return p

    def para(self, text, *, size=10.5, bold=False, italic=False, indent=None,
             space_after=4):
        p = self.doc.add_paragraph()
        p.paragraph_format.space_after = Pt(space_after)
        if indent is not None:
            p.paragraph_format.left_indent = Cm(indent)
        run = p.add_run(text)
        _style_run(run, size=size, bold=bold)
        run.font.italic = italic
        return p

    def bullet(self, text, size=10.5):
        p = self.doc.add_paragraph(style="List Bullet")
        p.paragraph_format.space_after = Pt(2)
        run = p.add_run(text)
        _style_run(run, size=size)
        return p

    def numbered(self, text, size=10.5):
        p = self.doc.add_paragraph(style="List Number")
        p.paragraph_format.space_after = Pt(2)
        run = p.add_run(text)
        _style_run(run, size=size)
        return p

    def quote(self, text):
        p = self.doc.add_paragraph(style="Quote")
        run = p.add_run(text)
        _style_run(run, size=10.5)
        return p

    def note(self, label, text):
        """Callout-style note paragraph."""
        p = self.doc.add_paragraph()
        p.paragraph_format.space_after = Pt(6)
        p.paragraph_format.left_indent = Cm(0.4)
        r1 = p.add_run(f"【{label}】")
        _style_run(r1, size=10.5, bold=True, color=RGBColor(0xC0, 0x00, 0x00))
        r2 = p.add_run(text)
        _style_run(r2, size=10.5)
        return p

    def table_caption(self, text):
        self.table_no += 1
        p = self.doc.add_paragraph()
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        p.paragraph_format.space_before = Pt(6)
        p.paragraph_format.space_after = Pt(2)
        run = p.add_run(f"表 {self.table_no}　{text}")
        _style_run(run, size=9, bold=True, east_asia=CJK_FONT_HEAD,
                   ascii_font=ASCII_FONT_HEAD)
        return self.table_no

    def table(self, headers, rows, caption=None, widths=None):
        if caption:
            self.table_caption(caption)
        t = self.doc.add_table(rows=1, cols=len(headers))
        t.style = "Table Grid"
        t.alignment = WD_TABLE_ALIGNMENT.CENTER
        for i, h in enumerate(headers):
            cell = t.rows[0].cells[i]
            cell.text = ""
            run = cell.paragraphs[0].add_run(h)
            _style_run(run, size=9, bold=True, east_asia=CJK_FONT_HEAD,
                       ascii_font=ASCII_FONT_HEAD)
            cell.paragraphs[0].alignment = WD_ALIGN_PARAGRAPH.CENTER
            _shade_cell(cell, HEADER_FILL)
        for row in rows:
            cells = t.add_row().cells
            for i, value in enumerate(row):
                cells[i].text = ""
                run = cells[i].paragraphs[0].add_run(str(value))
                _style_run(run, size=9)
        if widths:
            for row in t.rows:
                for i, w in enumerate(widths):
                    row.cells[i].width = Cm(w)
        _repeat_header_row(t)
        self.doc.add_paragraph().paragraph_format.space_after = Pt(2)
        return t

    def kv_table(self, pairs, caption=None):
        return self.table(["项目", "内容"], pairs, caption=caption,
                          widths=[4.0, 12.0])


# ---------------------------------------------------------------------------
# Document setup
# ---------------------------------------------------------------------------

def _setup_styles(doc):
    normal = doc.styles["Normal"]
    _set_style_fonts(normal, ASCII_FONT_BODY, CJK_FONT_BODY, size=10.5)
    normal.paragraph_format.line_spacing = 1.25
    normal.paragraph_format.space_after = Pt(4)
    for name in ("Heading 1", "Heading 2", "Heading 3"):
        try:
            _set_style_fonts(doc.styles[name], ASCII_FONT_HEAD, CJK_FONT_HEAD)
        except KeyError:
            pass


def _setup_page(doc):
    section = doc.sections[0]
    section.page_width = Cm(21.0)
    section.page_height = Cm(29.7)
    section.left_margin = Cm(2.4)
    section.right_margin = Cm(2.4)
    section.top_margin = Cm(2.4)
    section.bottom_margin = Cm(2.2)
    section.different_first_page_header_footer = True

    # header (not on cover page)
    hp = section.header.paragraphs[0]
    hp.alignment = WD_ALIGN_PARAGRAPH.RIGHT
    run = hp.add_run(f"抗体人源化评估标准指南（USAN / WHO INN）　{DOC_VERSION}")
    _style_run(run, size=8, color=RGBColor(0x80, 0x80, 0x80))

    # footer with page numbers (not on cover page)
    fp = section.footer.paragraphs[0]
    fp.alignment = WD_ALIGN_PARAGRAPH.CENTER
    _style_run(fp.add_run("第 "), size=8, color=RGBColor(0x80, 0x80, 0x80))
    _add_field(fp, " PAGE ")
    _style_run(fp.add_run(" 页 / 共 "), size=8, color=RGBColor(0x80, 0x80, 0x80))
    _add_field(fp, " NUMPAGES ")
    _style_run(fp.add_run(" 页"), size=8, color=RGBColor(0x80, 0x80, 0x80))


def _setup_core_properties(doc, date_str):
    props = doc.core_properties
    props.title = "抗体人源化评估标准指南（USAN / WHO INN）"
    props.subject = "抗体人源化程度评估、USAN/INN 命名解读与报告撰写规范"
    props.author = "humanization pipeline"
    props.last_modified_by = "humanization pipeline"
    props.category = "standard guide"
    props.keywords = "USAN; WHO INN; antibody humanization; IMGT; Kabat; CDR grafting"
    props.comments = (
        f"Document ID {DOC_ID}, version {DOC_VERSION}; generated by "
        "tools/usan/generate_usan_report.py"
    )
    props.revision = 2
    try:
        stamp = datetime.datetime.strptime(date_str, "%Y-%m-%d")
        props.created = stamp
        props.modified = stamp
    except ValueError:
        pass


# ---------------------------------------------------------------------------
# Cover / front matter
# ---------------------------------------------------------------------------

def _cover(b, date_str):
    doc = b.doc
    for _ in range(3):
        doc.add_paragraph()
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    _style_run(p.add_run("抗体人源化评估标准指南"), size=26, bold=True,
               east_asia=CJK_FONT_HEAD, ascii_font=ASCII_FONT_HEAD)
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    _style_run(p.add_run("基于 USAN / WHO INN 标准与 IMGT 编号体系"), size=14,
               color=RGBColor(0x66, 0x66, 0x66))
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    _style_run(p.add_run(f"标准编号 {DOC_ID}　版本 {DOC_VERSION}　{date_str}"),
               size=11, color=RGBColor(0x66, 0x66, 0x66))
    doc.add_paragraph()

    b.kv_table([
        ("文档编号", DOC_ID),
        ("版本", DOC_VERSION),
        ("生成日期", date_str),
        ("适用范围", "Fab（VH+VL）与 VHH 抗体的人源化程度评估、命名解读与结果报告"),
        ("编制工具", "tools/usan/generate_usan_report.py（本管线自带生成器）"),
        ("主要依据", "WHO INN Working Doc. 22-542 (2022)；Lefranc IMGT (2003)；"
                     "Kabat (1991)；本管线 scripts/humanize/"),
        ("免责声明", "本指南描述学术与行业通行惯例，供内部评估与论文撰写参考；"
                     "官方命名以 WHO / AMA 审定结果为准，监管要求以现行法规文本为准"),
    ], caption="文档控制信息")

    b.table_caption("修订历史")
    b.table(["版本", "日期", "修订说明"], [
        ("v1.0", "2026-09-11", "初版"),
        (DOC_VERSION, date_str,
         "按 WHO INN 22-542 (2022) 更正新方案四词干（-tug/-bart/-ment/-mig）；"
         "更正 IMGT CDR 边界（VH/VL 相同，CDR1=27-38）；将同源性阈值统一标注为学术估计；"
         "修正全人源实例（pembrolizumab/atezolizumab/natalizumab 实为人源化）；"
         "软化监管表述；新增术语表、目录、页码、本管线输出映射与附录"),
    ], widths=[2.0, 2.6, 11.4])

    doc.add_page_break()


def _toc(b):
    b.heading("目录", 1)
    p = b.doc.add_paragraph()
    _add_field(p, ' TOC \\o "1-3" \\h \\z \\u ')
    b.note("提示", "在 Microsoft Word 中打开后，右键点击目录 → “更新域”，"
                   "即可生成 / 刷新带页码的目录。")
    b.doc.add_page_break()


# ---------------------------------------------------------------------------
# Chapters
# ---------------------------------------------------------------------------

def _chapter_intro(b):
    b.heading("1　引言与适用范围", 1)

    b.heading("1.1　目的", 2)
    b.para("本指南规定本项目对抗体人源化程度的评估口径、USAN / WHO INN 抗体命名的解读方式，"
           "以及在论文、海报和报告中的规范表述，确保不同人、不同批次的结果可比、可复核。")
    b.para("本指南适用于：", space_after=2)
    b.bullet("Fab（VH+VL 配对）与 VHH 单域抗体的人源化程度评估；")
    b.bullet("CDR 移植 + 回复突变流程产生的 V0-V3 变体的报告；")
    b.bullet("对外材料（论文、海报、专利实施例、客户报告）中的命名与同源性表述。")

    b.heading("1.2　术语与缩写", 2)
    b.table(["术语", "说明"], GLOSSARY, caption="术语与缩写表", widths=[3.2, 12.8])

    b.heading("1.3　局限性与免责声明", 2)
    b.para("请在使用任何数值结论前阅读本节：")
    b.bullet("官方 USAN / INN 命名依据申办方申报的制备技术与分子格式（来源与工程化方式），"
             "并由 WHO / AMA 审定；不存在官方认可的“同源性百分比阈值”。")
    b.bullet("本指南第 3.5 节的 FR 同源性阈值（95 / 85 / 70 %）是学术文献中常用的近似估计，"
             "仅用于在只有序列、没有工艺信息时给出“估计类别”，不能作为命名依据。")
    b.bullet("FR 同源性仅覆盖 FR1+FR2+FR3：不含 CDR3（供体移植而来）、不含 FR4（J 基因编码）、"
             "不含恒定区，因此不能代表“整个抗体的鼠源序列比例”。")
    b.bullet("同源性是序列层面的单一维度指标：不反映三维结构、结合亲和力、稳定性与可开发性，"
             "评估结论需结合结构验证与实验数据。")
    b.bullet("不同编号体系（IMGT / Kabat / Chothia）的同源性数值不可直接互比；跨文献比较时必须注明编号体系。")


def _chapter_naming(b):
    b.heading("2　USAN / WHO INN 抗体命名体系", 1)

    b.heading("2.1　USAN 与 INN 的关系", 2)
    b.para("USAN（United States Adopted Name）由美国医学协会（AMA）药品命名委员会"
           "（USAN Council，历史上曾与美国药典 USP 共同管理）审定，作为美国药品的非专有名。")
    b.para("INN（International Nonproprietary Names）由 WHO 的 INN 计划主持，通过定期的 "
           "INN 例会与 Proposed / Recommended INN 名单发布。WHO 官方说明：各国非专有名"
           "（BAN、DCF、JAN、USAN 等）现已与 INN 基本一致。")
    b.para("实务要点：", space_after=2)
    b.bullet("命名由申办方向相应机构申请、审定后使用；商品名不得与 INN 常用词干混淆。")
    b.bullet("美国上市药品的标签使用非专有名（通常即 USAN / INN）；本指南不将任何命名机构表述为“审批机构”。")
    b.bullet("同一分子的名称一经确定即固定，后续序列微调（如回复突变）不改变其命名。")

    b.heading("2.2　传统来源后缀（旧方案，仍被广泛解读）", 2)
    b.para("旧方案以 V 区来源定义后缀，是文献中解读已获批药物名称的主要依据：")
    b.table(["后缀", "含义", "制备来源", "代表药物"], SUFFIX_TABLE,
            caption="传统抗体命名后缀（旧方案）", widths=[2.2, 3.4, 5.2, 5.2])
    b.note("注意", "后缀描述的是分子的来源类别（由申办方在申报中确定），"
                   "不是按序列同源性百分比判定；表中含义用于解读历史药物名称。")

    b.heading("2.3　2022 年 INN 新方案（按分子格式分类）", 2)
    b.para("2021 年 10 月第 73 届 INN 例会决定：停用统一词干 -mab，改用四个按“分子格式”"
           "划分的新词干；首批采用新方案的名称见 Proposed INN List 126"
           "（WHO INN Working Document 22-542, 2022）。")
    b.table(["词干", "类别", "定义要点", "对应分子形式"], NEW_STEMS_TABLE,
            caption="2022 年 INN 抗体新词干（WHO 22-542）", widths=[1.8, 3.2, 7.0, 4.0])
    b.para("新方案在词干之前使用“靶点中缀（infix）”标示作用类别，例如 -ta-（肿瘤）、"
           "-ci-（心血管）、-ne-（神经）等：")
    b.table(["中缀", "target class", "中文", "中缀", "target class", "中文"],
            [tuple(list(INFIX_TABLE[i]) + (list(INFIX_TABLE[i + 1])
                                           if i + 1 < len(INFIX_TABLE)
                                           else ["—", "—", "—"]))
             for i in range(0, len(INFIX_TABLE), 2)],
            caption="2022 年新方案靶点中缀（节选）",
            widths=[1.6, 4.0, 2.2, 1.6, 4.0, 2.2])
    b.para("方案补充说明：", space_after=2)
    b.bullet("抗体药物偶联物（ADC）同样适用本方案，不额外添加后缀（名称第二词即表明为偶联物）。")
    b.bullet("若融合蛋白的另一结构域并非免疫球蛋白来源（如细胞因子融合），归 -fusp 等其他方案，"
             "不适用抗体词干。")
    b.bullet("同一 CDR 组介导的交叉反应（如 bimekizumab）不属于 -mig。")

    b.heading("2.4　命名与序列同源性的关系（关键澄清）", 2)
    b.para("INN / USAN 命名依据分子的来源与格式类别（由申办方在申请中申报的制备技术与序列构成），"
           "WHO / AMA 不按序列同源性百分比分配后缀。因此：")
    b.numbered("拿到一个已获批药物名称，按 2.2 / 2.3 节解读其来源或格式类别；")
    b.numbered("只有序列、没有工艺信息时，可用 3.5 节的阈值给出“估计类别”，并明确标注为学术估计；")
    b.numbered("对外表述时始终区分“官方命名”与“学术估计”，避免写成“根据 FR 同源性被命名为 …”。")


def _chapter_eval(b):
    b.heading("3　人源化程度评估方法", 1)

    b.heading("3.1　重链与轻链分别评估", 2)
    b.para("VH 与 VL 必须分别评估、分别报告：", space_after=2)
    b.bullet("重链与轻链可能匹配不同的 germline 基因；")
    b.bullet("两条链的人源化程度通常不同，合并平均会掩盖短板；")
    b.bullet("文献与行业报告普遍分链给出 FR / CDR 同源性。")

    b.heading("3.2　CDR 分区定义", 2)
    b.para("推荐使用 IMGT 编号（Lefranc 2003）。IMGT 的关键性质：VH 与 VL 使用完全相同的区域边界；"
           "Kabat 编号下 VH 与 VL 的边界不同，引用时必须注明体系。")
    b.table(["CDR", "IMGT 位置", "Kabat 位置（VH / VL）", "本管线建议"], CDR_TABLE,
            caption="CDR 边界：IMGT 与 Kabat 对照", widths=[2.6, 3.4, 4.6, 5.4])
    b.note("更正说明", "旧版指南曾将 IMGT CDR-H1 写为 26-35、将轻链 IMGT 边界写成 Kabat 数值；"
                       "正确值为 IMGT CDR1 = 27-38（VH/VL 相同）。本文与 "
                       "scripts/humanize/imgt_numbering.py 的 IMGT_REGIONS 一致。")
    b.quote("“CDR 定义采用 IMGT 编号系统 (Lefranc, 2003)。”")

    b.heading("3.3　评估区域", 2)
    b.para("本管线的评估区域为 FR1 + FR2 + FR3（IMGT 1-26 / 39-55 / 66-104），排除 CDR3 与 FR4：")
    b.bullet("CDR3 由供体移植而来，且长度多变，计入会稀释框架同源性；")
    b.bullet("FR4 由人 J 基因编码（如 IGHJ1-6），天然高度保守，对区分人源化来源贡献很小；")
    b.bullet("IMGT 评估集合 = FR1 + CDR1 + FR2 + CDR2 + FR3（参见 imgt_numbering.py 的 IMGT_EVAL_REGIONS）。")
    b.note("定位", "上述区域划分是学术与行业报告的通行惯例（便于跨文献比较），"
                   "并非 USAN / INN 的官方规定——命名不基于该百分比。")

    b.heading("3.4　计算定义", 2)
    b.quote("FR 同源性 = FR1+FR2+FR3 中相同残基之和 / FR1+FR2+FR3 中可比较残基总数")
    b.quote("CDR 同源性 = CDR1+CDR2 中相同残基之和 / CDR1+CDR2 中可比较残基总数（不含 CDR3）")
    b.para("实现细节（与 scripts/humanize/imgt_numbering.py 的 compare_imgt_posmaps_direct 一致）：")
    b.bullet("仅对查询序列与 germline 同时存在的位置计数；单侧残基（query-only / germline-only）"
             "作为覆盖度诊断单独报告，不进入分母；")
    b.bullet("germline 为缺口（-）的位置跳过；")
    b.bullet("IMGT 编号下 VH / VL 区域边界相同，故 FR / CDR 统计口径两链一致。")

    b.heading("3.5　学术近似阈值（非官方）", 2)
    b.para("仅当无法获知制备工艺时，可用下表按 FR 同源性给出“估计类别”。"
           "对外表述必须附带估计字样：")
    b.table(["FR 同源性", "估计类别", "对应旧后缀", "典型来源"], THRESHOLD_TABLE,
            caption="学术近似阈值（依据 FR 同源性，IMGT）", widths=[3.4, 3.6, 3.2, 5.8])
    b.note("口径说明", "tools/humanness_score 按上表四档输出“Academic: likely …”标签；"
                       "humanization_result.json 的 usan_naming_class 采用简化的 85 / 70 两档"
                       "（≥0.85 Humanized、≥0.70 Chimeric、其余 Murine）。"
                       "两者均为学术估计，不是官方命名。")
    b.para("示例（规范写法）：")
    b.quote("“该抗体采用 CDR 移植技术人源化；基于 IMGT 编号，重链 FR 同源性 85.2%、"
            "轻链 FR 同源性 92.3%（均为与最匹配人 germline 的比较）。按学术惯例，该序列"
            "与人源化类别一致；具体命名以制备技术申报及 WHO / AMA 审定为准。”")


def _chapter_numbering(b):
    b.heading("4　编号体系：IMGT 与 Kabat", 1)

    b.heading("4.1　历史背景", 2)
    b.table(["体系", "创建时间", "创建者", "初始目的"], [
        ("Kabat", "1979 / 1991", "Elvin Kabat", "序列比对与结构研究"),
        ("IMGT", "1989 / 2003", "Marie-Paule Lefranc", "免疫遗传学数据库标准化"),
    ], caption="两种主要编号体系", widths=[2.4, 3.0, 4.4, 6.2])

    b.heading("4.2　核心差异", 2)
    b.bullet("Kabat：CDR 边界基于序列可变性；VH 与 VL 边界不同；插入用字母编码（如 H52A、H100B）。")
    b.bullet("IMGT：CDR 边界基于结构；VH 与 VL 使用相同位置定义；通过顶部缺口处理长度变化，"
             "CDR3 可延伸至 117 之后。")
    b.table(["区域", "IMGT 位置", "区域", "IMGT 位置", "区域", "IMGT 位置"],
            IMGT_REGION_TABLE, caption="IMGT 区域边界速查（VH = VL）",
            widths=[2.2, 2.8, 2.2, 2.8, 2.2, 2.8])

    b.heading("4.3　评估与工程化如何选择体系", 2)
    b.table(["用途", "建议体系", "理由"], [
        ("人源化程度评估 / 跨文献比较", "IMGT", "编号统一、数据库与多数文献采用，便于复核"),
        ("CDR 移植 / 回复突变设计", "Kabat（或与实验记录一致的体系）", "结构文献与历史工程数据以 Kabat 为主"),
        ("亲和力成熟位点标注", "与既有实验记录一致", "保持突变体与历史数据可对照"),
    ], caption="编号体系选择建议", widths=[5.2, 5.0, 5.8])
    b.note("重要", "同一份结果中不要混用体系；跨体系引用时给出换算表。"
                   "本管线 IMGT 模式（--numbering imgt）要求 ANARCI；Kabat 模式为便携默认。")
    b.para("本管线工作流：评估用 IMGT（标准化、便于比较），工程化标注用与实验记录一致的编号；"
           "两者互补，报告中必须注明各自体系。")


def _chapter_writing(b):
    b.heading("5　论文与报告撰写规范", 1)

    b.heading("5.1　命名与同源性的表述", 2)
    b.para("推荐句式（中文）：", space_after=2)
    b.quote("摘要：“该抗体采用 CDR 移植技术人源化。基于 IMGT 编号系统 (Lefranc, 2003)，"
            "与最匹配人 germline 比较，重链 FR 同源性为 85.2%，轻链 FR 同源性为 92.3%。"
            "命名后缀依据制备技术（WHO / INN 标准）确定。”")
    b.quote("方法：“人源化程度评估采用 IMGT 编号系统 (Lefranc, 2003)。FR 同源性定义为 "
            "FR1+FR2+FR3（排除 CDR3 与 FR4）中与人 germline 相同残基占可比较残基的比例；"
            "仅序列、无工艺信息时给出的类别为学术估计（≥95% fully-human-like，"
            "≥85% humanized-like，≥70% chimeric-like）。命名依据 WHO / INN 制备技术与格式标准。”")
    b.para("推荐句式（English）：", space_after=2)
    b.quote("\"The antibody was humanized by CDR grafting. Framework region (FR1+FR2+FR3, "
            "CDR3 and FR4 excluded) identity to the best-matching human germline was 85.2% "
            "(VH) and 92.3% (VL) using IMGT numbering (Lefranc, 2003). INN/USAN suffix "
            "assignment follows the WHO nomenclature scheme based on the manufacturing "
            "route and molecular format.\"")
    b.para("避免的写法：", space_after=2)
    b.bullet("“根据 FR 同源性 ≥85%，该抗体被命名为人源化抗体（-zumab）。”——命名不基于百分比。")
    b.bullet("“USAN / FDA 要求 FR 同源性达到 85%。”——不存在该官方阈值。")
    b.bullet("不注明编号体系就给出同源性百分比。")

    b.heading("5.2　结果表模板", 2)
    b.table(["链", "Germline", "FR 同源性 (IMGT)", "CDR 同源性 (IMGT)", "人源化程度（学术估计）"], [
        ("VH", "IGHV1-46*01", "85.2%", "62.5%", "humanized-like（-zumab 一致）"),
        ("VL", "IGKV1-39*01", "92.3%", "78.9%", "humanized-like（-zumab 一致）"),
    ], caption="人源化评估结果表示例", widths=[1.4, 3.4, 3.4, 3.4, 4.4])
    b.para("表注应包含：编号体系、评估区域（FR1+FR2+FR3；CDR1+CDR2）、"
           "germline 数据来源与版本、类别为学术估计的声明。")

    b.heading("5.3　补充评估维度（可选）", 2)
    b.table(["评估维度", "方法 / 来源", "说明"], [
        ("CVI 同源性", "Boehringer Ingelheim, JBC 2024；本管线 config.py 位置集",
         "canonical + vernier + interface 残基同源性，与表达量 / 亲和力保留相关"),
        ("Germline 使用频率", "治疗性抗体使用统计（Adimab / OAS 风格；本管线 germline_frequency 模块）",
         "已获批与临床阶段抗体中的 V 基因使用比例"),
        ("结构验证", "AF3 / 晶体结构；Kabsch 叠加后 CDR RMSD",
         "埋藏性与接触分析，指导回复突变分级"),
    ], caption="补充评估维度", widths=[3.0, 6.0, 7.0])


def _chapter_examples(b):
    b.heading("6　已获批抗体分类实例", 1)

    b.para("下表按旧方案后缀归类（-umab 全人源 / -zumab 人源化）。"
           "注意：部分公众熟知的药物常被误列为“全人源”，其后缀实际为人源化。")

    b.heading("6.1　全人源（后缀 -umab）", 2)
    b.table(["药物", "靶点", "适应症", "获批年份", "制备技术"], FULLY_HUMAN_EXAMPLES,
            caption="已批准的全人源抗体（-umab）", widths=[3.2, 2.4, 3.6, 2.2, 4.6])

    b.heading("6.2　人源化（后缀 -zumab）——易混淆对照", 2)
    b.table(["药物", "靶点", "适应症", "获批年份", "制备技术"], HUMANIZED_EXAMPLES,
            caption="已批准的人源化抗体（-zumab，含常被误列者）",
            widths=[3.2, 2.4, 3.6, 2.2, 4.6])
    b.note("易错点", "Natalizumab、Pembrolizumab、Atezolizumab 的后缀均为 -zumab，"
                     "属于人源化抗体，不应列入“全人源”表格。")

    b.heading("6.3　四类抗体对比", 2)
    b.table(["方面", "鼠源 (-omab)", "嵌合 (-ximab)", "人源化 (-zumab)", "全人源 (-umab)"], [
        ("V 区来源", "鼠", "鼠 V 区", "人框架 + 鼠 CDR", "人抗体库"),
        ("典型 FR 同源性", "< 70%", "70-85%", "85-95%", "≥ 95%"),
        ("典型 CDR 同源性", "低", "低-中", "中（供体 CDR）", "高"),
        ("免疫原性风险", "高", "中", "低-中", "低"),
        ("代表药物", "Muromonab", "Rituximab", "Trastuzumab", "Adalimumab"),
    ], caption="四类抗体特征对比（典型值，个体差异大）",
            widths=[3.0, 3.0, 3.0, 3.5, 3.5])
    b.note("说明", "FR / CDR 同源性为常见范围，仅用于定性对照；免疫原性还取决于具体残基、"
                   "暴露程度与个体差异，需实测或工具评估（如本管线 immunogenicity 模块）。")


def _chapter_mapping(b):
    b.heading("7　本管线输出与本指南的对应关系", 1)
    b.para("本指南是本仓库工具链的文字标准；下表把流水线产物映射到相应章节，便于复核。")
    b.table(["产物 / 字段", "含义", "对应章节"], JSON_MAPPING_TABLE,
            caption="本管线输出映射", widths=[5.4, 7.6, 3.0])
    b.para("标准命令：", space_after=2)
    b.bullet("同源性与命名估计（独立工具）：python3 tools/humanness_score/evaluate_humanness.py --input <fasta>")
    b.bullet("完整人源化流水线：python3 scripts/humanize/cli.py run --input <fasta> --outdir outputs")
    b.bullet("策略对比（Step 1）：python3 scripts/humanize/cli.py compare --input <fasta>")


def _chapter_refs(b):
    b.heading("8　参考文献", 1)
    for i, ref in enumerate(REFERENCES, 1):
        p = b.doc.add_paragraph()
        p.paragraph_format.space_after = Pt(2)
        p.paragraph_format.left_indent = Cm(0.75)
        p.paragraph_format.first_line_indent = Cm(-0.75)
        _style_run(p.add_run(f"[{i}]　{ref}"), size=9.5)


def _appendix(b):
    b.heading("附录 A　IMGT / Kabat 边界速查", 1)
    b.para("IMGT（VH 与 VL 相同）：FR1 1-26；CDR1 27-38；FR2 39-55；CDR2 56-65；"
           "FR3 66-104；CDR3 105-117（可延伸）；FR4 118-128。")
    b.para("Kabat VH：FR1 1-30；CDR1 31-35；FR2 36-49；CDR2 50-65；FR3 66-94；"
           "CDR3 95-102；FR4 103-113。")
    b.para("Kabat VL：FR1 1-23；CDR1 24-34；FR2 35-49；CDR2 50-56；FR3 57-88；"
           "CDR3 89-97；FR4 98-107。")
    b.para("实现位置：IMGT_REGIONS / KABAT_VH_REGIONS / KABAT_VL_REGIONS 定义于 "
           "scripts/humanize/imgt_numbering.py；评估集合 IMGT_EVAL_REGIONS = "
           "FR1+CDR1+FR2+CDR2+FR3（排除 CDR3 与 FR4）。")
    b.para("生成信息：本文件由 tools/usan/generate_usan_report.py 生成；内容如需修订，"
           "请修改该脚本后重新生成，避免直接编辑 .docx。")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def create_document(date_str: str = ""):
    if not date_str:
        date_str = datetime.date.today().isoformat()
    doc = Document()
    _setup_styles(doc)
    _setup_page(doc)
    _setup_core_properties(doc, date_str)
    b = GuideBuilder(doc)
    _cover(b, date_str)
    _toc(b)
    _chapter_intro(b)
    _chapter_naming(b)
    _chapter_eval(b)
    _chapter_numbering(b)
    _chapter_writing(b)
    _chapter_examples(b)
    _chapter_mapping(b)
    _chapter_refs(b)
    _appendix(b)
    return doc


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    parser = argparse.ArgumentParser(
        description="Generate the USAN / WHO INN antibody evaluation standard guide")
    parser.add_argument("--output", "-o",
                        default=os.path.join(here, "USAN标准指南.docx"),
                        help="output .docx path (default: alongside this script)")
    parser.add_argument("--date", default="",
                        help="pin the generation date (YYYY-MM-DD); default: today")
    args = parser.parse_args()

    out_dir = os.path.dirname(os.path.abspath(args.output))
    os.makedirs(out_dir, exist_ok=True)
    doc = create_document(args.date)
    doc.save(args.output)
    print(f"Word document generated: {args.output}")


if __name__ == "__main__":
    main()
