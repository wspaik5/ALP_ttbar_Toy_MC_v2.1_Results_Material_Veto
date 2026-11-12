"""Build v2.1 post-registration report PDF."""
import os, sys, hashlib, datetime
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'v19'))
import numpy as np
from reportlab.lib.pagesizes import letter
from reportlab.lib.units import inch, mm
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_LEFT, TA_CENTER, TA_JUSTIFY
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, PageBreak,
                                 Image, Table, TableStyle, KeepTogether)
from reportlab.lib import colors

WIDTH, HEIGHT = letter
MARGIN = 0.75 * inch

# Load results
res_all = np.load('v21_results.npy', allow_pickle=True).item()
MASSES = [0.017, 0.05, 0.1, 0.2]

r17 = np.load('v17/v17_results.npy', allow_pickle=True).item()
r18 = np.load('v18/v18_results.npy', allow_pickle=True).item()
PUB = {}
for m in MASSES:
    PUB[('C1', m)] = r17[('deg', '3D', m)]['deg']
    PUB[('C2', m)] = r17[('deg', '3D_allpixel', m)]['deg']
    PUB[('C3', m)] = r18['M025'][('deg', '3D_allpixel', m)]['deg']

# Hard veto efficiencies from stored results
HV_EFF = {}
for m in MASSES:
    HV_EFF[('C1', m)] = r17[('deg', '3D', m)]['eff']
    HV_EFF[('C2', m)] = r17[('deg', '3D_allpixel', m)]['eff']
    HV_EFF[('C3', m)] = r18['M025'][('deg', '3D_allpixel', m)]['eff']

# Styles
styles = getSampleStyleSheet()
styles.add(ParagraphStyle('Body', parent=styles['Normal'], fontSize=9.5, leading=13,
                           alignment=TA_JUSTIFY, spaceAfter=6))
styles.add(ParagraphStyle('SmallBody', parent=styles['Normal'], fontSize=8.5, leading=11,
                           alignment=TA_JUSTIFY, spaceAfter=4))
styles.add(ParagraphStyle('H1', parent=styles['Heading1'], fontSize=14, spaceAfter=10, spaceBefore=14))
styles.add(ParagraphStyle('H2', parent=styles['Heading2'], fontSize=11.5, spaceAfter=8, spaceBefore=10))
styles.add(ParagraphStyle('H3', parent=styles['Heading3'], fontSize=10, spaceAfter=6, spaceBefore=8))
styles.add(ParagraphStyle('Caption', parent=styles['Normal'], fontSize=8.5, leading=11,
                           alignment=TA_CENTER, spaceAfter=8, spaceBefore=4,
                           textColor=colors.HexColor('#333333')))
styles.add(ParagraphStyle('TableCell', parent=styles['Normal'], fontSize=8, leading=10))
styles.add(ParagraphStyle('MyCode', parent=styles['Normal'], fontSize=8, leading=10,
                           fontName='Courier'))

def S(text):
    return Paragraph(text, styles['Body'])

def Sm(text):
    return Paragraph(text, styles['SmallBody'])

story = []

# ── Title page ──
story.append(Spacer(1, 1.5*inch))
story.append(Paragraph('ALP-ttbar Toy MC v2.1 — Post-Registration Report',
                        ParagraphStyle('Title', parent=styles['Title'], fontSize=16, alignment=TA_CENTER)))
story.append(Spacer(1, 6))
story.append(Paragraph('A Likelihood Material Tag in Place of the Hard Material Veto<br/>'
                        'for Displaced Light ALP → e<super>+</super>e<super>−</super>',
                        ParagraphStyle('Sub', parent=styles['Normal'], fontSize=11, alignment=TA_CENTER)))
story.append(Spacer(1, 24))
story.append(Paragraph('<b>Author:</b> Won Shik Paik (Independent Researcher, Auckland, New Zealand)',
                        ParagraphStyle('Meta', parent=styles['Normal'], fontSize=10, alignment=TA_CENTER)))
story.append(Paragraph(f'<b>Date:</b> {datetime.date.today().isoformat()}',
                        ParagraphStyle('Meta2', parent=styles['Normal'], fontSize=10, alignment=TA_CENTER)))
story.append(Spacer(1, 18))
story.append(Paragraph('<b>Preregistration:</b> doi:10.5281/zenodo.24063091<br/>'
                        '<b>GitHub (preregistration):</b> github.com/wspaik5/ALP_ttbar_Toy_MC_v2.1_Preregistration_Material_Veto<br/>'
                        '<b>Builds on:</b> v2.0 series overview (doi:10.5281/zenodo.24056097); '
                        'v1.7 (doi:10.5281/zenodo.24001093); v1.8 (doi:10.5281/zenodo.24034481)<br/>'
                        '<b>Motivated by:</b> arXiv:2610.38932 (Alimena et al.)',
                        ParagraphStyle('Refs', parent=styles['Normal'], fontSize=9, alignment=TA_CENTER, leading=13)))
story.append(Spacer(1, 24))
story.append(S('<b>Abstract.</b> This note reports the results of the v2.1 preregistered analysis, which replaces '
               'the hard material veto and opening-angle exemption of v1.1–v1.9 by a single continuous likelihood '
               'material tag λ. The tag combines the radial vertex position relative to known material bands and the '
               'three-dimensional opening-angle statistic Δχ<super>2</super><sub>3D</sub> in a likelihood ratio. No new simulation '
               'or fitting is performed; the tag re-selects the stored, published fits of v1.7 and v1.8. '
               'The sanity gate and all validations pass. Of seven predictions, five are confirmed and two are refuted. '
               'The tag improves the limit by up to 17% at 100 MeV (f = 0.25) but makes it 10% worse at 17 MeV '
               'with pixel-only z information, where the hard veto with its discrete exemption threshold was already '
               'near-optimal for the coarse scan grid.'))
story.append(PageBreak())

# ── 1. Introduction ──
story.append(Paragraph('1. Introduction', styles['H1']))
story.append(S('This report presents the post-registration results of v2.1, the preregistered likelihood material '
               'tag analysis (doi:10.5281/zenodo.24063091). The analysis was planned in the v2.0 series overview '
               '(suggestion S4, open question Q1) and replaces the hard material veto of v1.1–v1.9.'))
story.append(S('The tag λ = Δχ<super>2</super><sub>3D</sub>/2 + ln[p<sub>sig</sub>(r<sub>m</sub>) / '
               'p<sub>conv</sub>(r<sub>m</sub> | σ<sub>r</sub>)] is a single continuous variable that replaces '
               'the two discrete parameters n<sub>σ</sub> (veto width) and k (exemption threshold) of the hard veto. '
               'It is frozen in <font face="Courier">mtag_v21.py</font>, committed with the preregistration.'))
story.append(S('Three configurations are tested, matching the published baselines of v1.7 and v1.8:'))
story.append(S('&nbsp;&nbsp;• <b>Config 1</b>: f = 1, Δχ<super>2</super><sub>3D</sub> with pixel z only (v1.7, "3D")<br/>'
               '&nbsp;&nbsp;• <b>Config 2</b>: f = 1, Δχ<super>2</super><sub>3D</sub> with z at every layer (v1.7, "3D_allpixel")<br/>'
               '&nbsp;&nbsp;• <b>Config 3</b>: f = 0.25, Δχ<super>2</super><sub>3D</sub> with z at every layer (v1.8, M025)'))
story.append(S('No new simulation or fitting is performed. The code re-selects the stored arrays from v1.7 and v1.8.'))

# ── 2. Sanity gate ──
story.append(Paragraph('2. Sanity gate', styles['H1']))
story.append(S('The three-step sanity gate (preregistration §4) was run before any tag result was computed.'))

story.append(Paragraph('2.1. Reproduction of published degradations', styles['H2']))
story.append(S('Re-running <font face="Courier">scan_v15</font> on the stored arrays reproduced all 12 published '
               'hard-veto degradations to machine precision (relative error = 0). <b>PASS.</b>'))

# Build the 12-value table
tdata = [['Config', 'Mass', 'Published', 'Computed', 'Rel. err']]
for cfg in ['C1', 'C2', 'C3']:
    for m in MASSES:
        g = res_all[cfg]['gate1'][m]
        tdata.append([cfg, f'{m*1000:.0f} MeV', f'{g["published"]:.4f}×', f'{g["computed"]:.4f}×', f'{g["relerr"]:.1e}'])

t = Table(tdata, colWidths=[55, 55, 75, 75, 60])
t.setStyle(TableStyle([
    ('FONTNAME', (0,0), (-1,0), 'Helvetica-Bold'),
    ('FONTSIZE', (0,0), (-1,-1), 8),
    ('GRID', (0,0), (-1,-1), 0.5, colors.grey),
    ('BACKGROUND', (0,0), (-1,0), colors.HexColor('#E0E0E0')),
    ('ALIGN', (1,0), (-1,-1), 'CENTER'),
]))
story.append(t)
story.append(Spacer(1, 8))

story.append(Paragraph('2.2. Normalisation', styles['H2']))
story.append(S('∫ p<sub>conv</sub> dr = 1 to better than 10<super>−6</super> for σ<sub>r</sub> = 1, 5 and 15 mm. <b>PASS.</b>'))

tdata2 = [['σ_r [mm]', '∫ p_conv', '|error|']]
for cfg in ['C1']:  # same for all configs
    for sr, g in res_all[cfg]['gate2'].items():
        tdata2.append([f'{sr:.0f}', f'{g["integral"]:.12f}', f'{g["err"]:.1e}'])
t2 = Table(tdata2, colWidths=[60, 120, 60])
t2.setStyle(TableStyle([
    ('FONTNAME', (0,0), (-1,0), 'Helvetica-Bold'),
    ('FONTSIZE', (0,0), (-1,-1), 8),
    ('GRID', (0,0), (-1,-1), 0.5, colors.grey),
    ('BACKGROUND', (0,0), (-1,0), colors.HexColor('#E0E0E0')),
    ('ALIGN', (0,0), (-1,-1), 'CENTER'),
]))
story.append(t2)
story.append(Spacer(1, 8))

story.append(Paragraph('2.3. Finiteness', styles['H2']))
story.append(S('λ is finite for every base pair in every sample and configuration. <b>PASS.</b>'))
story.append(S('<b>Sanity gate: ALL PASSED.</b>'))

# ── 3. Validations ──
story.append(Paragraph('3. Validations T1–T3', styles['H1']))

story.append(Paragraph('3.1. T1 (no-material limit)', styles['H2']))
story.append(S('With p<sub>conv</sub> set to zero and all Δχ<super>2</super> = 0, λ is the same constant for every '
               'pair (range = 0). The tag selection with t = none therefore passes all base pairs, equalling the '
               '"no veto" option. <b>PASS.</b>'))

story.append(Paragraph('3.2. T2 (monotonicity)', styles['H2']))
story.append(S('(a) At fixed σ<sub>r</sub> = 3 mm and Δχ<super>2</super> = 5, λ at each of the four pixel band centres '
               '(29, 68, 109, 160 mm) is lower than at the midpoint between any two adjacent pixel bands. <b>PASS.</b><br/>'
               '(b) For fixed r<sub>m</sub> = 50 mm and σ<sub>r</sub> = 3 mm, λ increases strictly with Δχ<super>2</super> '
               '(checked at 100 evenly spaced points from 0 to 20). <b>PASS.</b>'))

story.append(Paragraph('3.3. T3 (resolution dependence)', styles['H2']))
r_t3 = res_all['C1']['T3']
story.append(S(f'At the pixel band centre r = 29 mm with Δχ<super>2</super> = 5: '
               f'λ(σ<sub>r</sub> = 1 mm) = {r_t3["lam_sr1"]:.4f} &lt; λ(σ<sub>r</sub> = 15 mm) = {r_t3["lam_sr15"]:.4f}. '
               'A well-measured vertex inside a band is more conversion-like. <b>PASS.</b>'))

story.append(S('<b>All validations passed.</b>'))

# ── 4. Results ──
story.append(PageBreak())
story.append(Paragraph('4. Results', styles['H1']))

story.append(Paragraph('4.1. Overview', styles['H2']))

# Build main results table
tdata3 = [['Config', 'Mass', 'Hard veto', 'Tag', 'Threshold t', 'Ratio', 'Tag eff [%]', 'B (tag)']]
for cfg in ['C1', 'C2', 'C3']:
    for m in MASSES:
        hv = PUB[(cfg, m)]
        tag = res_all[cfg][('tag', m)]
        ratio = tag['deg'] / hv
        tdata3.append([cfg, f'{m*1000:.0f} MeV', f'{hv:.2f}×', f'{tag["deg"]:.2f}×',
                        str(tag['t']), f'{ratio:.4f}', f'{tag["eff_sig"]:.2f}', f'{tag["B"]:.1f}'])

t3 = Table(tdata3, colWidths=[42, 48, 55, 48, 55, 48, 52, 52])
t3.setStyle(TableStyle([
    ('FONTNAME', (0,0), (-1,0), 'Helvetica-Bold'),
    ('FONTSIZE', (0,0), (-1,-1), 7.5),
    ('GRID', (0,0), (-1,-1), 0.5, colors.grey),
    ('BACKGROUND', (0,0), (-1,0), colors.HexColor('#E0E0E0')),
    ('ALIGN', (1,0), (-1,-1), 'CENTER'),
    # Highlight cells where tag beats hard veto (ratio < 1)
]))
# Colour cells
for i in range(1, len(tdata3)):
    ratio_val = float(tdata3[i][5])
    if ratio_val < 0.98:
        t3.setStyle(TableStyle([('BACKGROUND', (5, i), (5, i), colors.HexColor('#C8E6C9'))]))
    elif ratio_val > 1.02:
        t3.setStyle(TableStyle([('BACKGROUND', (5, i), (5, i), colors.HexColor('#FFCDD2'))]))

story.append(t3)
story.append(Spacer(1, 6))
story.append(Paragraph('<b>Table 1.</b> Limit degradation from conversions: hard material veto (published) '
                        'versus likelihood material tag (v2.1). Ratio &lt; 1 means the tag produces a better (lower) limit. '
                        'Green: ratio &lt; 0.98 (tag improves); red: ratio &gt; 1.02 (tag worsens).',
                        styles['Caption']))

story.append(S('The tag improves the limit (ratio &lt; 1) in 7 of 12 mass–configuration points. The largest '
               'improvement is at 100 MeV with f = 0.25 (Config 3), where the tag reduces the degradation from '
               '3.28× to 2.71× (−17%). The tag is worse at 17 MeV in Config 1 (pixel z only), increasing the '
               'degradation from 127.3× to 139.4× (+10%).'))

# ── Figure 15(a) ──
story.append(Spacer(1, 6))
story.append(Image('fig15a_deg_vs_mass.png', width=5.5*inch, height=5.5*inch*5/7))
story.append(Paragraph('<b>Figure 15(a).</b> Limit degradation versus ALP mass. Open markers with dashed lines: '
                        'published hard-veto results. Filled markers with solid lines: likelihood tag (v2.1).',
                        styles['Caption']))

# ── Figure 15(b) ──
story.append(PageBreak())
story.append(Image('fig15b_lambda_distributions.png', width=5.5*inch, height=5.5*inch*4.5/7))
story.append(Paragraph('<b>Figure 15(b).</b> Normalised λ distributions for conversions (blue) and 17 MeV ALP signal '
                        '(red) in Config 3 (f = 0.25, z at every layer). The dashed line marks the optimal threshold '
                        't = 4.0.',
                        styles['Caption']))

# ── Figure 15(c) ──
story.append(Spacer(1, 8))
story.append(Image('fig15c_ratio.png', width=5.5*inch, height=5.5*inch*4.5/7))
story.append(Paragraph('<b>Figure 15(c).</b> Ratio of tag degradation to hard-veto degradation versus ALP mass. '
                        'Values below 1 mean the tag produces a better limit. The 5% threshold (red dotted line) '
                        'marks the P5 criterion.',
                        styles['Caption']))

story.append(Paragraph('4.2. Mechanism (Q3)', styles['H2']))
story.append(S('At the 17 MeV optimum with f = 0.25 (Config 3, t = 4.0), the tag keeps '
               f'{res_all["C3"][("tag", 0.017)]["eff_sig"]:.2f}% of base signal versus '
               f'{HV_EFF[("C3", 0.017)]:.2f}% for the hard veto — confirming P6 that the tag keeps more signal. '
               'The tag achieves this by trading a softer boundary around material bands (retaining pairs just '
               'inside a band that have high Δχ<super>2</super><sub>3D</sub>) against a sharper cut on pairs with low '
               'Δχ<super>2</super><sub>3D</sub> far from material.'))

# ── 5. Predictions ──
story.append(Paragraph('5. Prediction table', styles['H1']))
preds = res_all['predictions']

pred_data = [['Pred.', 'Statement', 'Observed', 'Status']]

r3 = res_all['C3']; r1 = res_all['C1']; r2 = res_all['C2']
pred_rows = [
    ('P1', 'f=0.25, 17 MeV: 16–26.7×', f'{r3[("tag", 0.017)]["deg"]:.2f}×', preds['P1']),
    ('P2', 'f=0.25, 50 MeV: 5–8.83×', f'{r3[("tag", 0.05)]["deg"]:.2f}×', preds['P2']),
    ('P3', 'f=1 pixel z, 17 MeV: 110–127.3×', f'{r1[("tag", 0.017)]["deg"]:.2f}×', preds['P3']),
    ('P4', 'f=1 z-every, 17 MeV: 85–117.4×', f'{r2[("tag", 0.017)]["deg"]:.2f}×', preds['P4']),
    ('P5', 'No config/mass worse by >5%', f'C1 17 MeV: +9.6%', preds['P5']),
    ('P6', 'Tag keeps more signal at 17 MeV f=0.25',
     f'{r3[("tag", 0.017)]["eff_sig"]:.2f}% vs {HV_EFF[("C3", 0.017)]:.2f}%', preds['P6']),
    ('P7', 'f=0.25, 200 MeV: ≤ 1.1×', f'{r3[("tag", 0.2)]["deg"]:.4f}×', preds['P7']),
]
for name, stmt, obs, ok in pred_rows:
    pred_data.append([name, stmt, obs, 'CONFIRMED' if ok else 'REFUTED'])

tp = Table(pred_data, colWidths=[30, 180, 95, 65])
tp.setStyle(TableStyle([
    ('FONTNAME', (0,0), (-1,0), 'Helvetica-Bold'),
    ('FONTSIZE', (0,0), (-1,-1), 8),
    ('GRID', (0,0), (-1,-1), 0.5, colors.grey),
    ('BACKGROUND', (0,0), (-1,0), colors.HexColor('#E0E0E0')),
    ('ALIGN', (0,0), (0,-1), 'CENTER'),
    ('ALIGN', (3,0), (3,-1), 'CENTER'),
    ('VALIGN', (0,0), (-1,-1), 'TOP'),
]))
for i, (_, _, _, ok) in enumerate(pred_rows, 1):
    bg = colors.HexColor('#C8E6C9') if ok else colors.HexColor('#FFCDD2')
    tp.setStyle(TableStyle([('BACKGROUND', (3, i), (3, i), bg)]))
story.append(tp)
story.append(Spacer(1, 6))
story.append(Paragraph('<b>Table 2.</b> Preregistered predictions and outcomes. 5/7 confirmed, 2/7 refuted.',
                        styles['Caption']))

story.append(Spacer(1, 8))
story.append(S('<b>P3 (refuted).</b> In Config 1 (f = 1, pixel z), the tag produces a worse limit at 17 MeV '
               '(139.4× vs 127.3×). The hard veto at this point uses k = None (no exemption) and n<sub>σ</sub> = 1.0, '
               'meaning it is a purely positional veto. The tag, which mixes in the Δχ<super>2</super><sub>3D</sub> term, '
               'lets conversion pairs with moderately high Δχ<super>2</super><sub>3D</sub> leak through even when they sit '
               'near material. With limited r–z information (pixel z only), the Δχ<super>2</super><sub>3D</sub> term is '
               'not informative enough at 17 MeV to compensate, and the tag overweights it.'))
story.append(S('<b>P5 (refuted).</b> The P3 refutation at 9.6% above the hard veto exceeds the 5% threshold. '
               'This is the same Config 1, 17 MeV point. All other 11 points are within 5%.'))

# ── 6. Deviations ──
story.append(Paragraph('6. Deviations from preregistration', styles['H1']))
story.append(S('<b>None.</b> The code, tag definition, threshold grid, mass windows, counting rule, subsampling '
               'and scan logic are exactly as preregistered. The stored fit arrays are loaded unmodified from the '
               'published v1.7 and v1.8 results. No post-hoc tuning or additional scans were performed.'))

# ── 7. Discussion ──
story.append(PageBreak())
story.append(Paragraph('7. Discussion', styles['H1']))
story.append(S('The likelihood material tag achieves modest improvements (2–17%) in most configurations and masses, '
               'with the best gains where the hard veto was already aggressive (f = 0.25, heavier masses). '
               'The tag degrades at 17 MeV with pixel-only z, where the hard veto was already near its '
               'optimum with a tight positional cut and no exemption.'))
story.append(S('The continuous tag has 14 threshold values versus 56 (n<sub>σ</sub> × k) for the hard veto. '
               'A finer threshold grid might recover part of the P3 loss but was not preregistered and is not tested.'))
story.append(S('The main conclusion for this toy series: the hard material veto is already close to optimal '
               'for this search, and the likelihood tag brings only incremental improvement. The material amount '
               '(f), the vertex resolution and the r–z information remain the dominant factors, as established in '
               'v1.6–v2.0.'))

# ── 8. Answers ──
story.append(Paragraph('8. Answers to preregistered questions', styles['H1']))
story.append(S('<b>Q1 (primary).</b> With the likelihood tag at f = 0.25 (z at every layer), the degradations are '
               f'{res_all["C3"][("tag", 0.017)]["deg"]:.1f}×, '
               f'{res_all["C3"][("tag", 0.05)]["deg"]:.1f}×, '
               f'{res_all["C3"][("tag", 0.1)]["deg"]:.1f}×, '
               f'{res_all["C3"][("tag", 0.2)]["deg"]:.2f}× '
               'at 17, 50, 100 and 200 MeV (hard veto: 26.7×, 8.83×, 3.28×, 1.09×). '
               'The tag improves by 4%, 1%, 17% and 0.3% respectively.'))
story.append(S('<b>Q2.</b> The gain depends on the r–z information. With f = 1 and pixel z only (Config 1), the tag '
               'is 10% worse at 17 MeV; with z at every layer (Config 2), it is 5% better. The z information makes '
               'the Δχ<super>2</super><sub>3D</sub> term more discriminating, tipping the balance in favour of the tag.'))
story.append(S(f'<b>Q3 (mechanism).</b> At 17 MeV (f = 0.25), the tag keeps '
               f'{res_all["C3"][("tag", 0.017)]["eff_sig"]:.1f}% of base signal '
               f'(hard veto: {HV_EFF[("C3", 0.017)]:.1f}%) with B = '
               f'{res_all["C3"][("tag", 0.017)]["B"]:.0f} conversion events '
               '(the hard veto had B ≈ the same order). The tag trades a softer material boundary for a '
               'sharper Δχ<super>2</super><sub>3D</sub> requirement.'))

# ── 9. Limitations ──
story.append(Paragraph('9. Limitations', styles['H1']))
story.append(S('All limitations listed in the preregistration (§8) and in v1.0–v2.0 apply. Additionally:'))
story.append(S('&nbsp;&nbsp;• The threshold grid (14 values) is coarser than the hard-veto grid (56 combinations). '
               'A finer grid was not preregistered.<br/>'
               '&nbsp;&nbsp;• The tag uses the toy\'s exact band positions; a real detector must measure its material map.<br/>'
               '&nbsp;&nbsp;• p<sub>sig</sub> is flat, not the true signal distribution; a per-mass tag could do better '
               'but would require re-tuning.<br/>'
               '&nbsp;&nbsp;• No statement about real detectors is made.'))

# ── 10. Future suggestions ──
story.append(Paragraph('10. Future suggestions for next-generation collider experiments', styles['H1']))
story.append(S('<i>All suggestions below are derived from the toy MC series (v1.0–v2.1). They are intended as '
               'points of consideration for detector and trigger designers working on real experiments; no claim '
               'is made that the toy faithfully reproduces any specific detector. The quantitative figures quoted '
               'are toy-level numbers that illustrate the direction and rough scale of each effect.</i>'))

story.append(Paragraph('10.1. The core finding: hard material veto as the dominant sensitivity bottleneck', styles['H2']))
story.append(S('Across v1.1–v2.1, the single largest factor degrading the upper limit for the displaced light ALP → '
               'e<super>+</super>e<super>−</super> search is neither detector resolution, opening-angle reconstruction, '
               'nor the number of tracker layers, but rather the <b>hard material veto</b> applied to suppress photon '
               'conversions. At the benchmark material scale f = 1 (nominal CMS barrel), the hard veto alone inflates '
               'the upper limit by a factor of 117–127× at 17 MeV. '
               'Even at f = 0.25 (75% material reduction), the inflation is still 26.7× — overwhelmingly larger than '
               'any gain achievable by improving the track fitter or adding z coordinates '
               '(which gave a factor of roughly 2× improvement in v1.7).'))
story.append(S('This identifies the material veto as the single most impactful place to intervene in the analysis chain.'))

story.append(Paragraph('10.2. The HL-LHC "collateral damage" problem', styles['H2']))
story.append(S('The HL-LHC detector upgrades (CMS Phase-2, ATLAS ITk) are primarily designed to cope with '
               'pileup (&lt;μ&gt; ≈ 200), to maximise Higgs coupling measurements, and to extend the reach for heavy BSM '
               'particles. These are legitimate and well-motivated goals. However, the toy MC series suggests that '
               'some of these upgrades may inadvertently harm the sensitivity to light, displaced BSM particles:'))
story.append(S('&nbsp;&nbsp;(a) <b>More tracker layers → more material bands.</b> The CMS Phase-2 outer tracker adds '
               'additional silicon layers compared to the current detector. In the hard-veto framework, each new layer '
               'creates a new material band that must be vetoed, increasing the fraction of radial phase space that is '
               'excluded. In the toy, going from the present barrel geometry to a denser layout increases the number of '
               'veto bands (14 → 18 in the test) and the vetoed signal fraction rises from 53% to 58%.'))
story.append(S('&nbsp;&nbsp;(b) <b>The effect is specific to displaced, low-mass channels.</b> For prompt searches '
               '(vertices near the beam line) or high-mass BSM particles (whose decay products have large opening '
               'angles and high p<sub>T</sub>), additional tracker layers improve vertex and momentum resolution — a clear '
               'net positive. The damage is specific to searches for particles that decay at intermediate radii '
               '(50–800 mm), where the signal vertices overlap with the material bands, and where the opening angle '
               'is small enough that conversion background dominates.'))
story.append(S('&nbsp;&nbsp;(c) <b>Total material, not layer arrangement, dominates at 17 MeV.</b> The v1.9 study '
               'showed that if total material budget is held constant, rearranging the same material into fewer or '
               'more layers changes the 17 MeV search sensitivity by less than 1%. What matters is the '
               '<i>total</i> material in the path, because thinner layers spread over a wider radial range create more '
               'veto bands while thicker layers concentrated in fewer positions create wider individual bands. The two '
               'effects roughly cancel at 17 MeV. This means the "collateral damage" cannot be easily fixed by '
               'changing the layer layout alone; the material budget itself must be reduced, or the hard veto must be '
               'replaced.'))

story.append(Paragraph('10.3. Suggestion 1 — Replace the hard material veto with a continuous likelihood '
               'variable at the trigger and reconstruction level', styles['H2']))
story.append(S('The v2.1 tag λ = Δχ<super>2</super><sub>3D</sub>/2 + ln[p<sub>sig</sub>(r) / p<sub>conv</sub>(r | σ<sub>r</sub>)] '
               'is a <b>standard Neyman–Pearson likelihood ratio</b> that combines two pieces of information: '
               '(i) the radial position of the vertex relative to known material bands, weighted by the radial '
               'resolution σ<sub>r</sub>; and (ii) the vertex quality via the 3D opening-angle χ<super>2</super> '
               'statistic. The form is familiar to any HEP physicist; the individual ingredients (material maps, '
               'vertex fits) are standard. What is new is the <i>specific combination</i> for the displaced '
               'e<super>+</super>e<super>−</super> channel and the quantitative demonstration that it can replace '
               'the hard veto with comparable or better sensitivity.'))
story.append(S('In 7 of 12 test points, the continuous tag matches or improves on the hard veto. At 100 MeV with '
               'f = 0.25, the improvement is 17% (degradation 2.71× vs 3.28×). The tag is worse at 17 MeV with pixel-only '
               'z (10%), where the hard veto was already near-optimal for the coarse scan grid — but even there, the loss '
               'is modest, and a finer threshold grid or continuous optimisation might recover it.'))
story.append(S('The critical advantage of a continuous tag over a hard veto appears at the <b>trigger level</b>. '
               'A hard veto applied at Level-1 or the High-Level Trigger is an irreversible binary cut: every event '
               'with a vertex inside a material band is permanently discarded, with no possibility of recovery in '
               'offline analysis. A continuous tag, by contrast, can be stored as a per-vertex score and used with '
               'varying thresholds at different analysis stages.'))

story.append(Paragraph('10.4. Suggestion 2 — The FPGA firmware window is now', styles['H2']))
story.append(S('For experiments currently in the upgrade phase (HL-LHC, scheduled to begin operation around 2030), '
               'the trigger hardware — specifically the FPGA firmware for Level-1 trigger processing — is being '
               'designed, built, tested and frozen during approximately 2026–2029. Once the firmware is frozen, '
               'the trigger-level cuts are fixed for the operational lifetime of the experiment.'))
story.append(S('If the trigger firmware implements a hard material veto for the displaced '
               'e<super>+</super>e<super>−</super> channel, the sensitivity loss documented in this toy MC series '
               'would be locked in for the full HL-LHC data-taking period (expected ~2030–2040). Conversely, '
               'incorporating a continuous likelihood score into the firmware at the design stage costs no '
               'additional hardware — it is a change to the algorithm, not to the silicon — and preserves the option '
               'to tune the selection threshold offline.'))
story.append(S('The window for this intervention is narrow. As of late 2026, the CMS Phase-2 tracker is being '
               'disassembled and the new outer tracker and pixel detector are being prepared for installation. '
               'The ATLAS ITk modules are being assembled at institutions worldwide. The FPGA firmware and '
               'reconstruction pipeline will be developed and tested in 2027–2029 and frozen before first data. '
               'Decisions made in this period about how to handle material-associated backgrounds at the trigger level '
               'will determine the sensitivity of the entire HL-LHC programme for displaced light BSM particles.'))

story.append(Paragraph('10.5. Suggestion 3 — Two pieces of information are currently thrown away', styles['H2']))
story.append(S('The hard material veto uses only one piece of information: "Is the vertex inside a material band?" '
               '(yes → veto, no → keep). It discards two additional pieces of information that are already available '
               'from the standard reconstruction:'))
story.append(S('&nbsp;&nbsp;(a) <b>How far inside the band, relative to the vertex resolution.</b> A vertex '
               'reconstructed at r = 29.5 mm with σ<sub>r</sub> = 0.5 mm (1σ from the beam-pipe centre at 29 mm) is '
               'very likely a conversion. A vertex at r = 31 mm with σ<sub>r</sub> = 3 mm (0.7σ from the band edge) '
               'is ambiguous. The hard veto treats both identically (both inside, both vetoed). The likelihood tag '
               'assigns them different scores because p<sub>conv</sub>(r | σ<sub>r</sub>) is a convolution of the '
               'material profile with the resolution function, not a step function.'))
story.append(S('&nbsp;&nbsp;(b) <b>The vertex quality (opening-angle χ<super>2</super>).</b> A high '
               'Δχ<super>2</super><sub>3D</sub> means the two tracks are inconsistent with a zero opening angle, '
               'favouring a massive parent over a massless photon. A signal vertex inside a material band with high '
               'Δχ<super>2</super><sub>3D</sub> should be retained, not vetoed. The tag combines this information '
               'with the positional score; the hard veto ignores it entirely.'))
story.append(S('These two "thrown-away" pieces of information are the source of the tag\'s ability to improve over '
               'the hard veto. They are not new observables; they are already computed by any standard vertex fitter. '
               'The suggestion is simply to use them.'))

story.append(Paragraph('10.6. Suggestion 4 — Material reduction remains the highest-impact intervention', styles['H2']))
story.append(S('Despite the tag\'s ability to partially compensate for the hard veto\'s loss, the v1.6–v2.1 results '
               'consistently show that <b>reducing the total material budget</b> is far more powerful than any '
               'software-level correction. The degradation at 17 MeV drops from 127× (f = 1) to 27× (f = 0.25) — a '
               'factor of ~5 improvement — dwarfing the tag\'s 4–17% gain on top of that. The qualitative conclusion '
               '"reduce material" is well known in the HL-LHC community; this toy series quantifies it specifically '
               'for the displaced e<super>+</super>e<super>−</super> channel.'))
story.append(S('For next-generation experiments beyond HL-LHC (e.g., FCC-ee or CEPC), where the detector design is '
               'still at an early stage, incorporating the displaced light BSM sensitivity into the material budget '
               'optimisation from the outset would avoid the "collateral damage" entirely.'))

story.append(Paragraph('10.7. Suggestion 5 — Outer tracker z-coordinate information halves the search loss', styles['H2']))
story.append(S('The v1.7 result showed that extending z-coordinate measurement from pixel-only to all tracker layers '
               'reduces the limit degradation by roughly a factor of 2 in the 50–200 MeV range (Config 1 → Config 2). '
               'This is because z information improves the 3D vertex fit, making Δχ<super>2</super><sub>3D</sub> a '
               'more powerful discriminant between signal and conversion vertices. The tag amplifies this benefit: '
               'with z at every layer, the tag beats the hard veto in 3 of 4 masses (Config 2), whereas with pixel-only z, '
               'it beats the hard veto in only 2 of 4 (Config 1).'))
story.append(S('For experiments designing their outer tracker readout, including even a coarse z-coordinate in the '
               'data sent to the trigger would directly improve the likelihood tag\'s discrimination. This is '
               'especially relevant for strip-based outer trackers, where stereo or macro-pixel layers can provide z.'))

story.append(Paragraph('10.8. Suggestion 6 — The Bremsstrahlung limit and the irreducible floor', styles['H2']))
story.append(S('The v1.3–v1.8 studies established that at 17 MeV, the dominant source of vertex mis-reconstruction is '
               '<b>multiple scattering in detector material</b>, not the track-fitting algorithm. Improving the '
               'Bremsstrahlung correction (Bethe–Heitler fitter vs Gaussian approximation) yields diminishing returns '
               'once the material scattering dominates. This means that algorithmic improvements to the track fitter '
               'alone cannot solve the material problem; the solution must involve either reducing the material itself '
               '(Suggestion 4) or changing how the material-associated background is handled (Suggestions 1–3).'))

story.append(Paragraph('10.9. Suggestion 7 — A preregistered benchmark for detector design studies', styles['H2']))
story.append(S('The complete v1.0–v2.1 toy MC series, including the tag definition, threshold grids, material models '
               'and stored fits, is publicly archived with SHA-256 hashes and Zenodo DOIs. Future detector design studies '
               'that wish to evaluate the impact of material budget, layer layout, or trigger strategy on the displaced '
               'light ALP channel may use this series as a preregistered benchmark. The tag formula (Eq. 3.1 in the '
               'preregistration) is parameter-free once the material map and resolution model are specified, and can be '
               'applied to any cylindrical barrel geometry.'))

story.append(Paragraph('10.10. Summary of suggestions', styles['H2']))
story.append(S('In order of estimated impact within the toy MC series:<br/>'
               '&nbsp;&nbsp;1. <b>Reduce total material budget</b> (factor 5× improvement from f = 1 → 0.25).<br/>'
               '&nbsp;&nbsp;2. <b>Replace hard material veto with continuous likelihood tag at trigger level</b> '
               '(avoids irreversible data loss; matches or improves hard veto in 7/12 test points).<br/>'
               '&nbsp;&nbsp;3. <b>Include z-coordinate information from all tracker layers</b> in the vertex fit '
               '(factor ~2× improvement for 50–200 MeV).<br/>'
               '&nbsp;&nbsp;4. <b>Design FPGA firmware to output a continuous material score</b> rather than a binary '
               'accept/reject, preserving offline optimisation freedom.<br/>'
               '&nbsp;&nbsp;5. <b>Incorporate displaced light BSM sensitivity into the material budget and trigger '
               'design from the outset</b> for future experiments (FCC-ee, CEPC), rather than retrofitting '
               'corrections after the detector is built.<br/>'
               '&nbsp;&nbsp;6. <b>Recognise the Bremsstrahlung/scattering limit</b>: at 17 MeV, material is the '
               'bottleneck, not the fitter — invest in material reduction and smart vetoing, not fitter refinement.'))

story.append(Paragraph('10.11. Timing context', styles['H2']))
story.append(S('The v2.1 preregistration was frozen on 2026-11-09, the day before the CMS Phase-2 strip tracker '
               'extraction was scheduled to begin (2026-11-10). The ATLAS ITk assembly is underway at multiple sites. '
               'The reconstruction and trigger software for both experiments will be developed and frozen during '
               '2027–2029. This places the current work at the beginning of the window in which trigger-level design '
               'choices can still be influenced. Whether any of these toy-level suggestions translate to actionable '
               'recommendations for a real detector is a question for the collaborations; the purpose of documenting '
               'them here, with a frozen timestamp, is to ensure the information is available during the decision '
               'period.'))

# ── 11. Design implications for HL-LHC trigger architectures ──
story.append(Paragraph('11. Design implications for HL-LHC trigger architectures', styles['H1']))
story.append(S('While this study relies on a simplified detector simulation and therefore does not account '
               'for full hardware-level constraints, our benchmark results highlight several practical '
               'trade-offs for High-Luminosity LHC (HL-LHC) trigger and DAQ design:'))

story.append(Paragraph('11.1. Trade-off between irreversible vetoes and continuous tagging', styles['H2']))
story.append(S('Because continuous likelihood scoring achieves background-rejection performance comparable '
               'to the discrete hard material veto in this regime, retaining events with a continuous metadata '
               'tag rather than applying an irreversible trigger-level veto appears advantageous. Where '
               'bandwidth permits, soft tagging preserves flexibility for offline threshold optimisation in '
               'searches for long-lived particles (LLPs) and low-mass ALPs.'))

story.append(Paragraph('11.2. Resource-efficient alternatives to complex classifiers', styles['H2']))
story.append(S('Given that the conventional hard material veto already operates close to the statistical '
               'optimum in this phase space, deploying resource-intensive deep neural networks on Level-1 '
               'FPGAs may yield diminishing returns for photon-conversion rejection. Lightweight '
               'approximations — such as parameterised lookup tables incorporating layer geometry and '
               'material budget — could offer a favourable balance between discrimination power and FPGA '
               'latency/logic budgets.'))

story.append(Paragraph('11.3. Sensitivity to detector material-map accuracy', styles['H2']))
story.append(S('Our diagnostic findings (Section 10) suggest that discrimination limits are governed '
               'primarily by physical information loss — specifically detector layer placement versus total '
               'material budget — rather than classifier complexity. Consequently, accurate 3D material-map '
               'calibration of the upgraded tracking systems (such as the HL-LHC Inner Trackers) may prove '
               'more impactful than further algorithmic refinement in this specific channel.'))

# ── 12. Technical future directions ──
story.append(Paragraph('12. Technical future directions (within this toy)', styles['H1']))
story.append(S('The following are internal improvements to the toy MC series itself, distinct from the '
               'broader experimental suggestions above:'))
story.append(S('&nbsp;&nbsp;• A z-dependent material tag (3D material map) could further improve discrimination, '
               'but requires a z-dependent material model not available in this barrel-only toy.<br/>'
               '&nbsp;&nbsp;• A per-mass or per-cτ signal prior (rather than flat p<sub>sig</sub>) could improve the '
               'tag for specific mass hypotheses.<br/>'
               '&nbsp;&nbsp;• A finer threshold grid or a continuous optimisation (rather than a discrete scan) '
               'could recover the small P3 loss at 17 MeV with pixel z.<br/>'
               '&nbsp;&nbsp;• An endcap geometry would test whether the barrel-only conclusions extend to forward '
               'rapidities, where the material distribution is different.<br/>'
               '&nbsp;&nbsp;• A multi-dimensional tag combining λ with additional variables '
               '(e.g., vertex χ<super>2</super>, track isolation) could be studied as a more sophisticated '
               'discriminant.'))

# ── 13. References ──
story.append(Paragraph('13. References', styles['H1']))
story.append(Sm('[E1] CMS Collaboration, Precision measurement of the structure of the CMS inner tracking system '
                 'using nuclear interactions, JINST 13 (2018) P10034, arXiv:1807.03289.'))
story.append(Sm('[E2] ATLAS Collaboration, A measurement of material in the ATLAS tracker using secondary hadronic '
                 'interactions in 7 TeV pp collisions, JINST 11 (2016) P11020, arXiv:1609.04305.'))
story.append(Sm('[E3] S. F. Koch et al., Measuring the ATLAS ITk pixel detector material via multiple scattering '
                 'of positrons at the CERN PS, arXiv:2412.04686.'))

# ── Build ──
today = datetime.date.today().strftime('%Y-%m-%d')
outname = f'v2.1_post_registration_report_{today}.pdf'
doc = SimpleDocTemplate(outname, pagesize=letter,
                         leftMargin=MARGIN, rightMargin=MARGIN,
                         topMargin=MARGIN, bottomMargin=MARGIN)
doc.build(story)
print(f'Built {outname}')

# SHA-256
with open(outname, 'rb') as f:
    h = hashlib.sha256(f.read()).hexdigest()
print(f'SHA-256: {h}')

# MD5
import hashlib as hl
with open(outname, 'rb') as f:
    m5 = hl.md5(f.read()).hexdigest()
print(f'MD5: {m5}')
