"""Patent Example Report Generation (专利实施例报告).

Generates Word document in patent example format with:
  - Humanization design description
  - Back-mutation design table
  - Variant sequences with SEQ ID NO
  - CDR analysis tables (Kabat and IMGT)
  - VH/VL cross-combination matrix

Requires python-docx.
"""

from __future__ import annotations

import os
from typing import Dict, List, Optional, Tuple

from .pipeline import ChainReport, RunResult


def _get_cdr_sequences_kabat(sequence: str, chain_type: str) -> Dict[str, str]:
    """Extract CDR sequences using Kabat numbering."""
    from .numbering import number_heavy, number_light
    
    if chain_type == 'H':
        numbered = number_heavy(sequence)
        kabat_cdrs = {
            'CDR1': ('H31', 'H35'),
            'CDR2': ('H50', 'H65'),
            'CDR3': ('H95', 'H102'),
        }
    else:
        numbered = number_light(sequence)
        kabat_cdrs = {
            'CDR1': ('L24', 'L34'),
            'CDR2': ('L50', 'L56'),
            'CDR3': ('L89', 'L97'),
        }
    
    posmap = numbered.posmap()
    
    def extract_cdr_seq(posmap, start_label, end_label):
        start_num = int(''.join(c for c in start_label if c.isdigit()))
        end_num = int(''.join(c for c in end_label if c.isdigit()))
        prefix = start_label[0]
        
        seq = ''
        for num in range(start_num, end_num + 1):
            label = f'{prefix}{num}'
            aa = posmap.get(label, '-')
            if aa != '-':
                seq += aa
            
            for letter in 'ABCDEFGHIJK':
                ins_label = f'{prefix}{num}{letter}'
                ins_aa = posmap.get(ins_label, None)
                if ins_aa:
                    seq += ins_aa
        
        return seq
    
    cdrs = {}
    for region, (start, end) in kabat_cdrs.items():
        cdrs[region] = extract_cdr_seq(posmap, start, end)
    
    return cdrs


def _get_cdr_sequences_imgt(sequence: str, chain_type: str) -> Dict[str, str]:
    """Extract CDR sequences using IMGT numbering via ANARCI.
    
    Uses ANARCI for proper IMGT numbering, then extracts CDR sequences
    based on IMGT region definitions (CDR1: 27-38, CDR2: 56-65, CDR3: 105-117).
    """
    from .imgt_numbering import IMGT_REGIONS
    
    # Try ANARCI first for accurate IMGT numbering
    try:
        from .anarci_adapter import is_anarci_available, number_with_anarci_imgt
        
        if is_anarci_available():
            numbered = number_with_anarci_imgt(sequence, chain_type)
            if numbered is not None:
                posmap = numbered.posmap()
                
                # Extract CDR sequences based on IMGT regions
                cdrs = {}
                prefix = chain_type
                for region, (start, end) in IMGT_REGIONS.items():
                    if region.startswith('CDR'):
                        seq = ''
                        for num in range(start, end + 1):
                            # Check for insertions (IMGT uses gaps, but ANARCI may have insertions)
                            label = f"{prefix}{num}"
                            aa = posmap.get(label, '-')
                            if aa != '-':
                                seq += aa
                            
                            # Check for insertions after this position
                            for letter in 'ABCDEFGHIJK':
                                ins_label = f"{prefix}{num}{letter}"
                                ins_aa = posmap.get(ins_label, None)
                                if ins_aa:
                                    seq += ins_aa
                        cdrs[region] = seq
                
                return cdrs
    except ImportError:
        pass
    except Exception as e:
        # Fallback to Kabat-based extraction
        pass
    
    # Fallback: Use Kabat numbering and convert
    from .numbering import number_heavy, number_light
    
    if chain_type == 'H':
        numbered = number_heavy(sequence)
    else:
        numbered = number_light(sequence)
    
    posmap = numbered.posmap()
    
    # Convert Kabat posmap to IMGT posmap
    from .imgt_numbering import kabat_posmap_to_imgt_posmap
    imgt_posmap = kabat_posmap_to_imgt_posmap(posmap, chain_type)
    
    # Extract CDR sequences
    cdrs = {}
    prefix = chain_type
    for region, (start, end) in IMGT_REGIONS.items():
        if region.startswith('CDR'):
            seq = ''
            for num in range(start, end + 1):
                label = f"{prefix}{num}"
                aa = imgt_posmap.get(label, '-')
                if aa != '-':
                    seq += aa
                
                for letter in 'ABCDEFGHIJK':
                    ins_label = f"{prefix}{num}{letter}"
                    ins_aa = imgt_posmap.get(ins_label, None)
                    if ins_aa:
                        seq += ins_aa
            cdrs[region] = seq
    
    return cdrs


def _get_backmutation_description(vh_backmuts: List[str], vl_backmuts: List[str],
                                   vh_candidates: List = None, vl_candidates: List = None) -> str:
    """Generate back-mutation description string.
    
    Args:
        vh_backmuts: List of VH back-mutation position strings (e.g., ["H67", "H29"])
        vl_backmuts: List of VL back-mutation position strings (e.g., ["L1"])
        vh_candidates: Optional list of BackMutationCandidate objects for VH
        vl_candidates: Optional list of BackMutationCandidate objects for VL
    """
    parts = []
    
    # Build position -> candidate mapping
    vh_map = {}
    if vh_candidates:
        for c in vh_candidates:
            vh_map[c.position] = c
    
    vl_map = {}
    if vl_candidates:
        for c in vl_candidates:
            vl_map[c.position] = c
    
    if vh_backmuts:
        vh_desc_parts = []
        for pos in vh_backmuts:
            if pos in vh_map:
                c = vh_map[pos]
                vh_desc_parts.append(f"{pos}{c.donor_aa}{c.human_aa}")
            else:
                vh_desc_parts.append(pos)
        vh_desc = ",".join(vh_desc_parts)
        parts.append(f"VH: {vh_desc}")
    
    if vl_backmuts:
        vl_desc_parts = []
        for pos in vl_backmuts:
            if pos in vl_map:
                c = vl_map[pos]
                vl_desc_parts.append(f"{pos}{c.donor_aa}{c.human_aa}")
            else:
                vl_desc_parts.append(pos)
        vl_desc = ",".join(vl_desc_parts)
        parts.append(f"VL: {vl_desc}")
    
    return "; ".join(parts) if parts else "None"


def build_patent_example_report(result: RunResult, out_path: str) -> str:
    """Build patent example Word document.
    
    Args:
        result: RunResult from humanization pipeline
        out_path: Output file path
    
    Returns:
        Path to generated document
    """
    from docx import Document
    from docx.enum.table import WD_TABLE_ALIGNMENT
    from docx.enum.text import WD_ALIGN_PARAGRAPH
    from docx.shared import Inches, Pt, RGBColor
    
    doc = Document()
    
    # ---- base styles ----
    style = doc.styles["Normal"]
    style.font.name = "Times New Roman"
    style.font.size = Pt(10.5)
    
    for hname, sz in [("Heading 1", 14), ("Heading 2", 12), ("Heading 3", 11)]:
        h = doc.styles[hname]
        h.font.name = "Times New Roman"
        h.font.size = Pt(sz)
    
    def para(text="", bold=False, italic=False, size=10.5, align=None, space_after=6):
        p = doc.add_paragraph()
        if align:
            p.alignment = align
        p.paragraph_format.space_after = Pt(space_after)
        r = p.add_run(text)
        r.bold, r.italic = bold, italic
        r.font.size = Pt(size)
        return p
    
    def seq_para(label, seq, size=9):
        p = doc.add_paragraph()
        p.paragraph_format.space_after = Pt(4)
        r = p.add_run(f"{label} ")
        r.bold = True
        r.font.size = Pt(size)
        r2 = p.add_run(seq)
        r2.font.name = "Courier New"
        r2.font.size = Pt(size)
        return p
    
    def table(headers, rows, widths=None, mono_cols=()):
        t = doc.add_table(rows=1 + len(rows), cols=len(headers))
        t.style = "Table Grid"
        t.alignment = WD_TABLE_ALIGNMENT.CENTER
        for j, h in enumerate(headers):
            c = t.cell(0, j)
            c.text = ""
            r = c.paragraphs[0].add_run(str(h))
            r.bold = True
            r.font.size = Pt(9)
        for i, row in enumerate(rows):
            for j, v in enumerate(row):
                c = t.cell(i + 1, j)
                c.text = ""
                r = c.paragraphs[0].add_run(str(v))
                r.font.size = Pt(8.5)
                if j in mono_cols:
                    r.font.name = "Courier New"
        if widths:
            for j, w in enumerate(widths):
                for row in t.rows:
                    row.cells[j].width = Inches(w)
        doc.add_paragraph()
        return t
    
    # ============================ TITLE ============================
    para("Example: Humanization Design of Antibodies", bold=True, size=14,
         align=WD_ALIGN_PARAGRAPH.CENTER, space_after=12)
    para("实施例：抗体的人源化设计", bold=True, size=12,
         align=WD_ALIGN_PARAGRAPH.CENTER, space_after=24)
    
    # ============================ DESIGN DESCRIPTION (EN) ============================
    para("Humanization design was performed through CDR grafting strategy using "
         "antibody humanization pipeline. By searching against the IMGT human antibody "
         "heavy and light chain variable region germline gene database, germline genes "
         "with high homology to the antibody variable regions were selected as candidate "
         "templates for the heavy and light chains respectively. The amino acid numbering "
         "and CDR regions were determined by the Kabat numbering scheme of the antibodies. "
         "The CDRs of the antibody were grafted into the framework regions (FR) of the "
         "corresponding human templates, forming variable region sequences in the order "
         "FR1-CDR1-FR2-CDR2-FR3-CDR3-FR4. The suitability of templates was comprehensively "
         "scored based on multiple factors including template gene homology, frequency of "
         "occurrence in clinical therapeutic antibodies and human antibodies, expression "
         "yield predictions, isoelectric point of post-grafting sequences, and high-risk "
         "PTM (post-translational modification) sites. High-scoring template sequences after "
         "grafting were preferentially selected as the basis for back mutation design.", 
         size=10.5, space_after=12)
    
    # ============================ DESIGN DESCRIPTION (CN) ============================
    para("使用抗体人源化流程的CDR移植策略进行人源化设计。通过比对IMGT人类抗体重轻链可变区"
         "种系基因数据库，分别挑选与抗体可变区同源性较高的重链和轻链可变区种系基因作为模板。"
         "使用Kabat编号系统确定抗体的氨基酸编号及CDR区域。将抗体的CDR分别移植到相应的模板"
         "的框架区（FR）中，形成次序为FR1-CDR1-FR2-CDR2-FR3-CDR3-FR4的可变区序列。通过"
         "计算模板基因的同源性、在临床及人体中的出现频率、历史表达量数据、移植后序列的等电点、"
         "PTM高风险位点等多重因素对模板的适用性进行综合评分，优选高分模板的移植后序列作为回复"
         "突变设计的基础。", size=10.5, space_after=12)
    
    # ============================ GERMLINE SELECTION ============================
    para("1. Germline Selection (种系基因选择)", bold=True, size=11, space_after=6)
    para("1. 种系基因选择", bold=True, size=11, space_after=6)
    
    # Collect germline info
    germline_info = []
    for rep in result.chains:
        c = rep.input_chain
        v = rep.germline.v_gene
        j = rep.germline.j_gene
        if v:
            germline_info.append({
                'chain': c.chain_type,
                'v_gene': v.gene_id,
                'j_gene': j.gene_id if j else 'N/A',
                'v_seq': v.sequence,
                'j_seq': j.sequence if j else '',
            })
    
    vh_germline = next((g for g in germline_info if g['chain'] == 'H'), None)
    vl_germline = next((g for g in germline_info if g['chain'] == 'L'), None)
    
    if vh_germline:
        para(f"The final selected humanization templates for the heavy chain were "
             f"{vh_germline['v_gene']} and {vh_germline['j_gene']}.", 
             size=10.5, space_after=6)
        para(f"最终选择的重链模板为{vh_germline['v_gene']}和{vh_germline['j_gene']}。",
             size=10.5, space_after=6)
    
    if vl_germline:
        para(f"The final selected humanization templates for the light chain were "
             f"{vl_germline['v_gene']} and {vl_germline['j_gene']}.",
             size=10.5, space_after=6)
        para(f"最终选择的轻链模板为{vl_germline['v_gene']}和{vl_germline['j_gene']}。",
             size=10.5, space_after=12)
    
    # ============================ BACK-MUTATION DESIGN ============================
    para("2. Back-mutation Design (回复突变设计)", bold=True, size=11, space_after=6)
    para("2. 回复突变设计", bold=True, size=11, space_after=6)
    
    para("Based on the predicted variable region structure, the impact of non-conservative "
         "amino acids in FR on CDR conformation was quantitatively assessed. Key amino acids "
         "with high scores were back mutated to the corresponding amino acids in the original "
         "antibody to maintain the original affinity. Different back mutation combinations "
         "resulted in multiple humanized variants. Selected back mutation designs are shown "
         "in Table 1.", size=10.5, space_after=6)
    
    para("基于预测的可变区结构，定量评估FR中的非保守氨基酸对CDR构象的影响，将评分较高的"
         "关键氨基酸回复突变为抗体对应的氨基酸，以维持原有的亲和力。不同的回复突变组合得到"
         "多个变体。具体回复突变设计见表1。", size=10.5, space_after=12)
    
    # Build back-mutation table
    para("Table 1: Back mutation design for humanized antibodies", bold=True, size=10, space_after=4)
    para("表1：人源化抗体回复突变设计", bold=True, size=10, space_after=6)
    
    # Collect variant info
    vh_variants = []
    vl_variants = []
    vh_backmut_candidates = []
    vl_backmut_candidates = []
    
    for rep in result.chains:
        c = rep.input_chain
        if c.chain_type == 'H':
            vh_variants = rep.variants
            vh_backmut_candidates = rep.backmut.candidates
        else:
            vl_variants = rep.variants
            vl_backmut_candidates = rep.backmut.candidates
    
    # Build position -> candidate mapping
    vh_map = {c.position: c for c in vh_backmut_candidates}
    vl_map = {c.position: c for c in vl_backmut_candidates}
    
    # Get V0 origin for comparison
    vh_v0_origin = vh_variants[0].graft.origin if vh_variants else {}
    vl_v0_origin = vl_variants[0].graft.origin if vl_variants else {}
    
    # Build table rows
    backmut_rows = []
    max_variants = max(len(vh_variants), len(vl_variants), 1)
    
    for i in range(max_variants):
        vl_desc = ""
        vh_desc = ""
        
        if i < len(vl_variants):
            v = vl_variants[i]
            # Use origin differences to match humanization_report logic
            diffs = [p for p in v.graft.origin
                     if v.graft.origin[p] != vl_v0_origin.get(p)]
            if diffs:
                bm_desc_parts = []
                for pos in sorted(diffs, key=lambda x: (x[0], int(''.join(c for c in x if c.isdigit())))):
                    if pos in vl_map:
                        c = vl_map[pos]
                        # 专利格式：位置+人源氨基酸+供体氨基酸 (e.g., L2IL)
                        # 供体插入位点：human_aa='-'，显示为 L6:E (供体插入E)
                        if c.human_aa and c.human_aa != '-':
                            bm_desc_parts.append(f"{pos}{c.human_aa}{c.donor_aa}")
                        elif c.donor_aa and c.donor_aa != '-':
                            # Donor insertion - show as pos:donor_aa
                            bm_desc_parts.append(f"{pos}:{c.donor_aa}")
                        else:
                            bm_desc_parts.append(pos)
                    else:
                        bm_desc_parts.append(pos)
                if bm_desc_parts:
                    bm_desc = ",".join(bm_desc_parts)
                    vl_desc = f"Graft({vl_germline['v_gene'] if vl_germline else 'VL'}) + {bm_desc}"
                else:
                    vl_desc = f"Graft({vl_germline['v_gene'] if vl_germline else 'VL'})"
            else:
                vl_desc = f"Graft({vl_germline['v_gene'] if vl_germline else 'VL'})"
        
        if i < len(vh_variants):
            v = vh_variants[i]
            # Use origin differences to match humanization_report logic
            diffs = [p for p in v.graft.origin
                     if v.graft.origin[p] != vh_v0_origin.get(p)]
            if diffs:
                bm_desc_parts = []
                for pos in sorted(diffs, key=lambda x: (x[0], int(''.join(c for c in x if c.isdigit())))):
                    if pos in vh_map:
                        c = vh_map[pos]
                        # 专利格式：位置+人源氨基酸+供体氨基酸 (e.g., H5VL)
                        # 供体插入位点：human_aa='-'，显示为 H6:E (供体插入E)
                        if c.human_aa and c.human_aa != '-':
                            bm_desc_parts.append(f"{pos}{c.human_aa}{c.donor_aa}")
                        elif c.donor_aa and c.donor_aa != '-':
                            # Donor insertion - show as pos:donor_aa
                            bm_desc_parts.append(f"{pos}:{c.donor_aa}")
                        else:
                            bm_desc_parts.append(pos)
                    else:
                        bm_desc_parts.append(pos)
                if bm_desc_parts:
                    bm_desc = ",".join(bm_desc_parts)
                    vh_desc = f"Graft({vh_germline['v_gene'] if vh_germline else 'VH'}) + {bm_desc}"
                else:
                    vh_desc = f"Graft({vh_germline['v_gene'] if vh_germline else 'VH'})"
            else:
                vh_desc = f"Graft({vh_germline['v_gene'] if vh_germline else 'VH'})"
        
        backmut_rows.append([f"L{i+1}" if i < len(vl_variants) else "", vl_desc,
                            f"H{i+1}" if i < len(vh_variants) else "", vh_desc])
    
    table(["VL", "", "VH", ""], backmut_rows, widths=[0.5, 2.5, 0.5, 3.0])
    
    para("Note: Graft represents grafting the antibody CDRs into the human germline "
         "template FR sequence; AxxxB indicates mutating A at position xxx to B, and so on "
         "(numbering follows positional index in grafted sequence).", 
         size=9, italic=True, space_after=6)
    
    para("注：Graft代表将抗体CDR植入人类模板FR区序列; AxxxB表示将第xxx位A突变成B，"
         "其它依此类推。回复突变氨基酸的编号为自然顺序编号。", 
         size=9, italic=True, space_after=12)
    
    # ============================ VARIANT SEQUENCES ============================
    para("3. Variant Sequences (变体序列)", bold=True, size=11, space_after=6)
    para("3. 变体序列", bold=True, size=11, space_after=6)
    
    para("The specific sequences of humanized antibody variable regions are as follows:",
         size=10.5, space_after=6)
    para("人源化抗体可变区具体序列如下：", size=10.5, space_after=12)
    
    # VL sequences
    para("Light chain variable region sequences:", bold=True, size=10, space_after=4)
    para("轻链可变区序列：", bold=True, size=10, space_after=6)
    
    seq_id_no = 1
    for i, v in enumerate(vl_variants):
        para(f"VL{i+1} amino acid sequence as shown in SEQ ID NO: {seq_id_no}",
             size=10, space_after=2)
        para(f"VL{i+1} 氨基酸序列如SEQ ID NO：{seq_id_no}所示：",
             size=10, space_after=4)
        seq_para("Sequence", v.sequence)
        seq_id_no += 1
    
    doc.add_paragraph()
    
    # VH sequences
    para("Heavy chain variable region sequences:", bold=True, size=10, space_after=4)
    para("重链可变区序列：", bold=True, size=10, space_after=6)
    
    for i, v in enumerate(vh_variants):
        para(f"VH{i+1} amino acid sequence as shown in SEQ ID NO: {seq_id_no}",
             size=10, space_after=2)
        para(f"VH{i+1} 氨基酸序列如SEQ ID NO：{seq_id_no}所示：",
             size=10, space_after=4)
        seq_para("Sequence", v.sequence)
        seq_id_no += 1
    
    doc.add_paragraph()
    
    # Template sequences
    para("Template sequences:", bold=True, size=10, space_after=4)
    para("模板序列：", bold=True, size=10, space_after=6)
    
    if vl_germline:
        para(f"Light chain V gene template {vl_germline['v_gene']} sequence as shown in "
             f"SEQ ID NO: {seq_id_no}", size=10, space_after=2)
        para(f"轻链模板{vl_germline['v_gene']}氨基酸序列如SEQ ID NO：{seq_id_no}所示：",
             size=10, space_after=4)
        seq_para("Sequence", vl_germline['v_seq'])
        seq_id_no += 1
        
        if vl_germline['j_seq']:
            para(f"Light chain FR4 template sequence as shown in SEQ ID NO: {seq_id_no}",
                 size=10, space_after=2)
            para(f"轻链FR4模板氨基酸序列如SEQ ID NO：{seq_id_no}所示：",
                 size=10, space_after=4)
            seq_para("Sequence", vl_germline['j_seq'])
            seq_id_no += 1
    
    doc.add_paragraph()
    
    if vh_germline:
        para(f"Heavy chain V gene template {vh_germline['v_gene']} sequence as shown in "
             f"SEQ ID NO: {seq_id_no}", size=10, space_after=2)
        para(f"重链模板{vh_germline['v_gene']}氨基酸序列如SEQ ID NO：{seq_id_no}所示：",
             size=10, space_after=4)
        seq_para("Sequence", vh_germline['v_seq'])
        seq_id_no += 1
        
        if vh_germline['j_seq']:
            para(f"Heavy chain FR4 template sequence as shown in SEQ ID NO: {seq_id_no}",
                 size=10, space_after=2)
            para(f"重链FR4模板氨基酸序列如SEQ ID NO：{seq_id_no}所示：",
                 size=10, space_after=4)
            seq_para("Sequence", vh_germline['j_seq'])
            seq_id_no += 1
    
    doc.add_page_break()
    
    # ============================ CROSS COMBINATION ============================
    para("4. Cross-combination of Variants (变体交叉组合)", bold=True, size=11, space_after=6)
    para("4. 变体交叉组合", bold=True, size=11, space_after=6)
    
    para("This invention selects different light chain and heavy chain sequences from "
         "the above back mutation designs of humanized antibody light chain and heavy "
         "chain variable regions for cross-combination, ultimately obtaining multiple "
         "humanized antibodies.", size=10.5, space_after=6)
    
    para("本发明分别从上述人源化抗体轻链和重链可变区的回复突变设计中，选择不同的轻链"
         "和重链序列进行交叉组合，最终获得多种人源化抗体。", size=10.5, space_after=12)
    
    para("Table 2: Antibody variable region corresponding amino acid sequences", 
         bold=True, size=10, space_after=4)
    para("表2：抗体可变区对应氨基酸序列", bold=True, size=10, space_after=6)
    
    # Build combination matrix
    combo_headers = ["Fv"] + [f"VH{i+1}" for i in range(len(vh_variants))]
    combo_rows = []
    for i, vl_v in enumerate(vl_variants):
        row = [f"VL{i+1}"]
        for j, vh_v in enumerate(vh_variants):
            row.append(f"Ab-{i*len(vh_variants)+j+1}")
        combo_rows.append(row)
    
    table(combo_headers, combo_rows)
    
    doc.add_page_break()
    
    # ============================ CDR ANALYSIS (KABAT) ============================
    para("5. CDR Analysis (CDR 分析)", bold=True, size=11, space_after=6)
    para("5. CDR 分析", bold=True, size=11, space_after=12)
    
    para("Based on the Kabat numbering scheme, the CDRs of the above humanized antibody "
         "VH and VL sequences are shown in Table 3.", size=10.5, space_after=6)
    para("根据Kabat编号系统，上述人源化抗体VH和VL序列分析结果如表3所示。",
         size=10.5, space_after=12)
    
    para("Table 3: CDRs of humanized antibody VH and VL sequences in Kabat",
         bold=True, size=10, space_after=4)
    para("表3：人源化抗体VH和VL序列的Kabat分析结果",
         bold=True, size=10, space_after=6)
    
    # Kabat CDR table
    kabat_headers = ["Variable", "CDR1", "CDR2", "CDR3"]
    kabat_rows = []
    
    for i, v in enumerate(vl_variants):
        cdrs = _get_cdr_sequences_kabat(v.sequence, 'L')
        kabat_rows.append([f"VL{i+1}", cdrs.get('CDR1', ''), cdrs.get('CDR2', ''), cdrs.get('CDR3', '')])
    
    for i, v in enumerate(vh_variants):
        cdrs = _get_cdr_sequences_kabat(v.sequence, 'H')
        kabat_rows.append([f"VH{i+1}", cdrs.get('CDR1', ''), cdrs.get('CDR2', ''), cdrs.get('CDR3', '')])
    
    table(kabat_headers, kabat_rows, mono_cols=(1, 2, 3))
    
    doc.add_page_break()
    
    # ============================ CDR ANALYSIS (IMGT) ============================
    para("Based on the IMGT numbering scheme, the CDRs of the above humanized antibody "
         "VH and VL sequences are shown in Table 4.", size=10.5, space_after=6)
    para("根据IMGT编号系统，上述人源化抗体VH和VL序列分析结果如表4所示。",
         size=10.5, space_after=12)
    
    para("Table 4: CDRs of humanized antibody VH and VL sequences in IMGT",
         bold=True, size=10, space_after=4)
    para("表4：人源化抗体VH和VL序列的IMGT分析结果",
         bold=True, size=10, space_after=6)
    
    # IMGT CDR table
    imgt_headers = ["Variable", "CDR1", "CDR2", "CDR3"]
    imgt_rows = []
    
    for i, v in enumerate(vl_variants):
        cdrs = _get_cdr_sequences_imgt(v.sequence, 'L')
        imgt_rows.append([f"VL{i+1}", cdrs.get('CDR1', ''), cdrs.get('CDR2', ''), cdrs.get('CDR3', '')])
    
    for i, v in enumerate(vh_variants):
        cdrs = _get_cdr_sequences_imgt(v.sequence, 'H')
        imgt_rows.append([f"VH{i+1}", cdrs.get('CDR1', ''), cdrs.get('CDR2', ''), cdrs.get('CDR3', '')])
    
    table(imgt_headers, imgt_rows, mono_cols=(1, 2, 3))
    
    # Save document
    doc.save(out_path)
    return out_path
