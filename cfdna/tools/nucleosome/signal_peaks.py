#!/usr/bin/env python3
"""
signal_peaks: Position/peak calling for continuous normalized coverage signals.

Reimplementation of core DANPOS algorithms adapted for per-base-pair
normalized coverage data in the range [-1, 1]. Replaces Poisson-based
statistics with z-score and variance-ratio approaches appropriate for
continuous signals. Removes rpy2 dependency entirely.

Designed for methylation difference, normalized coverage, or any
continuous per-bp signal in bigWig or wiggle format.

Usage:
    python signal_peaks.py dpos signal.bw [options]
    python signal_peaks.py dpos signal.wig [options]

Author: Adapted from DANPOS (Kaifu Chen et al.) for continuous signal data.
"""

import os
import sys
import argparse
import numpy as np
from time import time
from copy import deepcopy
from math import sqrt, log
from scipy import stats
import pyBigWig


# ---------------------------------------------------------------------------
# Wig class: Load, store, and manipulate per-chromosome signal arrays
# ---------------------------------------------------------------------------

class Wig:
    """
    Store genome-wide signal as {chrom: numpy_array}.

    Supports bigWig (.bw, .bigwig, .bigWig) and wiggle (.wig, .wig.gz) formats.
    Step=1 is the default for per-base-pair data; larger steps are supported
    for binned data.
    """

    BIGWIG_EXTENSIONS = ('.bw', '.bigwig', '.bigWig')

    def __init__(self, file="", step=0):
        self.data = {}   # {chrom_name: numpy_array}
        self.step = step
        if file:
            self.load(file)

    # ---- I/O ----

    def load(self, file):
        """Auto-detect bigWig or wiggle format by extension and load."""
        if any(file.endswith(ext) for ext in self.BIGWIG_EXTENSIONS):
            self._load_bigwig(file)
        else:
            self._load_wiggle(file)

    def _load_wiggle(self, file):
        """Detect fixed vs variable step wiggle and dispatch."""
        opener = open
        if file.endswith('.gz'):
            import gzip
            opener = gzip.open

        with opener(file, 'rt') as fi:
            for line in fi:
                line = line.strip()
                if not line or line.startswith('track'):
                    continue
                if line.startswith('fixedStep'):
                    self._load_fixed(file)
                    return
                elif line.startswith('variableStep'):
                    self._load_variable(file)
                    return
                else:
                    raise ValueError(f"Unrecognized wiggle format: {line[:40]}")

    def _load_bigwig(self, file):
        """
        Load signal from a bigWig file using pyBigWig.

        If self.step is set, the signal is binned at that resolution.
        If self.step == 0 (auto), step is set to 1 for per-base-pair
        resolution.
        """
        bw = pyBigWig.open(file)
        chroms = bw.chroms()

        if self.step < 1:
            self.step = 1

        step = self.step

        for chrom, chrom_len in chroms.items():
            n_bins = chrom_len // step
            if n_bins < 1:
                continue

            if step == 1:
                # Per-bp: use values() for exact signal
                vals = bw.values(chrom, 0, chrom_len)
                arr = np.array(vals, dtype=np.float64)
                arr = np.nan_to_num(arr, nan=0.0)
            else:
                # Binned: use stats() for mean over each bin
                arr = np.array(
                    bw.stats(chrom, 0, n_bins * step, nBins=n_bins, type="mean"),
                    dtype=np.float64
                )
                arr = np.nan_to_num(arr, nan=0.0)

            self.data[chrom] = arr

        bw.close()

    def _load_fixed(self, file):
        """Load fixed-step wiggle file."""
        opener = open
        if file.endswith('.gz'):
            import gzip
            opener = lambda f: gzip.open(f, 'rt')

        # Auto-detect step if not set
        if self.step < 1:
            with opener(file) as fi:
                for line in fi:
                    if line.startswith('fixedStep'):
                        for term in line.split()[1:]:
                            k, v = term.split('=')
                            if k == 'step':
                                self.step = int(v)
                                break
                    if self.step >= 1:
                        break
            if self.step < 1:
                self.step = 1

        chrom = None
        values = []
        start = 0
        instep = self.step

        with opener(file) as fi:
            for line in fi:
                line = line.strip()
                if not line or line.startswith('track'):
                    continue

                if line.startswith('fixedStep'):
                    # Flush previous chromosome block
                    if chrom is not None and values:
                        self._store_block(chrom, start, values, instep)
                        values = []

                    # Parse header
                    for term in line.split()[1:]:
                        k, v = term.split('=')
                        if k == 'chrom':
                            chrom = v
                        elif k == 'start':
                            start = int(v) - 1  # wiggle is 1-based
                        elif k == 'step':
                            instep = int(v)
                else:
                    try:
                        values.append(float(line.split()[0]))
                    except ValueError:
                        if line.startswith('nan'):
                            values.append(0.0)
                        else:
                            print(f"WARNING: skipping line: {line}", file=sys.stderr)

        # Flush last block
        if chrom is not None and values:
            self._store_block(chrom, start, values, instep)

    def _store_block(self, chrom, start, values, instep):
        """Store a block of values into the chromosome array."""
        arr = np.array(values, dtype=np.float64)

        # Resample if input step differs from target step
        if instep != self.step:
            if instep < self.step:
                # Downsample by averaging
                ratio = self.step // instep
                trim = len(arr) - len(arr) % ratio
                arr = arr[:trim].reshape(-1, ratio).mean(axis=1)
            else:
                # Upsample by repeating
                ratio = instep // self.step
                arr = np.repeat(arr, ratio)

        offset = start // self.step
        needed = offset + len(arr)

        if chrom not in self.data:
            self.data[chrom] = np.zeros(needed, dtype=np.float64)

        if self.data[chrom].size < needed:
            self.data[chrom] = np.resize(self.data[chrom], needed)

        self.data[chrom][offset:offset + len(arr)] = arr

    def _load_variable(self, file):
        """Load variable-step wiggle file."""
        if self.step < 1:
            self.step = 1

        opener = open
        if file.endswith('.gz'):
            import gzip
            opener = lambda f: gzip.open(f, 'rt')

        chrom = None
        span = 1

        with opener(file) as fi:
            for line in fi:
                line = line.strip()
                if not line or line.startswith('track'):
                    continue

                if line.startswith('variableStep'):
                    for term in line.split()[1:]:
                        k, v = term.split('=')
                        if k == 'chrom':
                            chrom = v
                            if chrom not in self.data:
                                self.data[chrom] = np.zeros(1, dtype=np.float64)
                        elif k == 'span':
                            span = int(v)
                else:
                    col = line.split()
                    pos = int(col[0]) - 1  # 1-based to 0-based
                    val = float(col[1])
                    s = pos // self.step
                    e = (pos + span) // self.step
                    if e >= self.data[chrom].size:
                        new_size = e + 10000
                        old = self.data[chrom]
                        self.data[chrom] = np.zeros(new_size, dtype=np.float64)
                        self.data[chrom][:old.size] = old
                    self.data[chrom][s:e] = val

        # Trim trailing zeros
        for chrom in self.data:
            last_nonzero = np.max(np.nonzero(self.data[chrom]), initial=0)
            self.data[chrom] = self.data[chrom][:last_nonzero + 1]

    def save(self, file, format='fixed'):
        """Save to wiggle format file."""
        step = self.step
        chroms = sorted(self.data.keys())
        with open(file, 'w') as f:
            for chrom in chroms:
                f.write(f"fixedStep chrom={chrom} start=1 step={step} span={step}\n")
                for val in self.data[chrom]:
                    f.write(f"{val}\n")
        print(f"Saved: {file}")

    # ---- Statistics ----

    def mean(self):
        total, count = 0.0, 0
        for chrom in self.data:
            total += self.data[chrom].sum()
            count += self.data[chrom].size
        return total / max(count, 1)

    def std(self):
        m = self.mean()
        ss, count = 0.0, 0
        for chrom in self.data:
            diff = self.data[chrom] - m
            ss += np.sum(diff * diff)
            count += self.data[chrom].size
        return sqrt(ss / max(count, 1))

    def gsize(self):
        return sum(self.data[c].size * self.step for c in self.data)

    def chr_size(self, chrom):
        if chrom in self.data:
            return self.data[chrom].size * self.step
        return 0

    # ---- Transforms ----

    def smooth(self, width=100):
        """Smooth signal with a uniform sliding window."""
        if width <= self.step:
            return
        w = max(1, width // self.step)
        kernel = np.ones(w) / w
        for chrom in self.data:
            self.data[chrom] = np.convolve(self.data[chrom], kernel, mode='same')

    def abs_signal(self):
        """Return a new Wig with absolute values (for bidirectional peak calling)."""
        out = Wig(step=self.step)
        for chrom in self.data:
            out.data[chrom] = np.abs(self.data[chrom])
        return out


# ---------------------------------------------------------------------------
# Summits class: Store and manipulate called summits
# ---------------------------------------------------------------------------

class Summits:
    """
    Store summit positions and values per chromosome.

    data[chrom] = {
        'p': numpy_array of summit positions (bp),
        'v': numpy_array of summit signal values,
        's': numpy_array of fuzziness scores (optional),
        'ppos': numpy_array of positioning p-values (optional),
        'zscore': numpy_array of z-scores (optional),
    }
    """

    def __init__(self):
        self.data = {}

    def fetch_values(self, wig):
        """Retrieve signal value at each summit position from a Wig."""
        for chrom in self.data:
            ps = self.data[chrom]['p']
            vals = np.zeros(ps.size, dtype=np.float64)
            if chrom in wig.data:
                for i in range(ps.size):
                    idx = ps[i] // wig.step
                    if idx < wig.data[chrom].size:
                        vals[i] = wig.data[chrom][int(idx)]
            self.data[chrom]['v'] = vals

    def compute_zscores(self, wig):
        """Compute z-score for each summit relative to genome-wide distribution."""
        m = wig.mean()
        s = wig.std()
        if s == 0:
            s = 1e-10
        for chrom in self.data:
            self.data[chrom]['zscore'] = (self.data[chrom]['v'] - m) / s

    def compute_fuzziness(self, wig, rd=100):
        """
        Compute fuzziness (positioning) score for each summit.

        Uses variance-ratio (F-test) comparing the observed read distribution
        around a summit to a uniform background distribution. Adapted from
        DANPOS log10fuztest but using scipy instead of rpy2.

        Parameters:
            wig: Wig instance with the signal
            rd: half-size of flanking region for fuzziness calculation (bp)
        """
        step = wig.step
        # Background variance: variance of positions under uniform distribution
        positions = np.arange(-rd, rd + step, step, dtype=np.float64)
        bg_var = np.var(positions)
        bg_n = len(positions)

        for chrom in self.data:
            ps = self.data[chrom]['p']
            n = ps.size
            fuz_scores = np.zeros(n, dtype=np.float64)
            fuz_pvals = np.zeros(n, dtype=np.float64)

            if chrom not in wig.data:
                self.data[chrom]['s'] = fuz_scores
                self.data[chrom]['ppos'] = fuz_pvals
                continue

            w = wig.data[chrom]
            max_idx = w.size - 1

            for i in range(n):
                center = ps[i] // step
                lo = int(center - rd // step)
                hi = int(center + rd // step + 1)

                if lo < 0 or hi > max_idx:
                    fuz_scores[i] = sqrt(bg_var)
                    fuz_pvals[i] = 0.0
                    continue

                # Weighted variance: weight each position's squared distance
                # by the signal value (shifted to be non-negative)
                region = w[lo:hi]
                weights = np.abs(region)  # use absolute signal as weights
                total_w = weights.sum()

                if total_w < 1e-10:
                    fuz_scores[i] = sqrt(bg_var)
                    fuz_pvals[i] = 0.0
                    continue

                pos_offsets = np.arange(lo - center, hi - center, dtype=np.float64) * step
                weighted_var = np.sum(weights * pos_offsets**2) / total_w
                obs_n = int(total_w + 0.5)
                obs_n = max(obs_n, 2)

                fuz_scores[i] = sqrt(weighted_var)

                if weighted_var >= bg_var:
                    fuz_pvals[i] = 0.0
                else:
                    f_stat = weighted_var / bg_var
                    # F-test: smaller variance in numerator
                    pval = stats.f.logcdf(f_stat, obs_n, bg_n)
                    fuz_pvals[i] = pval / np.log(10)  # convert to log10

            self.data[chrom]['s'] = fuz_scores
            self.data[chrom]['ppos'] = fuz_pvals

    def merge(self, wig, distance=110):
        """
        Merge neighboring summits closer than `distance` bp.

        When two summits are within distance, keep the one with the
        larger absolute value. If equal, take the midpoint.
        """
        step = wig.step
        for chrom in self.data:
            ps = self.data[chrom]['p']
            vs = self.data[chrom]['v']

            if ps.size < 2:
                continue

            merge_count = 1
            while merge_count > 0:
                merge_count = 0
                keep = []

                i = 0
                while i < ps.size:
                    if i == ps.size - 1:
                        keep.append(i)
                        break

                    if ps[i + 1] - ps[i] >= distance:
                        keep.append(i)
                        i += 1
                    else:
                        merge_count += 1
                        # Keep the summit with larger absolute value
                        if abs(vs[i]) > abs(vs[i + 1]):
                            vs[i + 1] = vs[i]
                            ps[i + 1] = ps[i]
                        elif abs(vs[i]) == abs(vs[i + 1]):
                            mid = (ps[i] + ps[i + 1]) // 2
                            ps[i + 1] = mid
                            idx = int(mid // step)
                            if chrom in wig.data and idx < wig.data[chrom].size:
                                vs[i + 1] = wig.data[chrom][idx]
                        # else: vs[i+1] already larger, skip i
                        i += 1

                if merge_count > 0:
                    keep_idx = np.array(keep, dtype=int)
                    # Also include the last element if it wasn't already
                    if len(keep) == 0 or keep[-1] != ps.size - 1:
                        keep_idx = np.append(keep_idx, ps.size - 1)
                    ps = ps[keep_idx]
                    vs = vs[keep_idx]

            self.data[chrom]['p'] = ps
            self.data[chrom]['v'] = vs

    def fillgap(self, wig, height=0.05, distance=110):
        """
        Insert a summit between two neighbors if gap size is appropriate.

        For continuous signals, `height` is the minimum absolute signal
        value required for a filled summit (default 0.05).
        """
        step = wig.step
        min_gap = distance * 2.5
        max_gap = distance * 3.5

        for chrom in self.data:
            ps = self.data[chrom]['p']
            vs = self.data[chrom]['v']

            if ps.size < 2:
                continue

            inserts_pos = []
            inserts_val = []
            inserts_idx = []

            for i in range(ps.size - 1):
                gap = ps[i + 1] - ps[i]
                if min_gap < gap < max_gap:
                    mid_bp = (ps[i] + ps[i + 1]) // 2
                    mid_idx = int(mid_bp // step)
                    if chrom in wig.data and mid_idx < wig.data[chrom].size:
                        val = wig.data[chrom][mid_idx]
                        if abs(val) >= height:
                            inserts_idx.append(i + 1)
                            inserts_pos.append(mid_bp)
                            inserts_val.append(val)

            if inserts_pos:
                new_ps = np.insert(ps, inserts_idx, inserts_pos)
                new_vs = np.insert(vs, inserts_idx, inserts_val)
                self.data[chrom]['p'] = new_ps
                self.data[chrom]['v'] = new_vs


# ---------------------------------------------------------------------------
# Summit calling: sliding-window approach adapted for continuous signals
# ---------------------------------------------------------------------------

def call_summits(wig, width=40, height=0.1, mode='abs'):
    """
    Call signal summits (local maxima) using a sliding window.

    Parameters:
        wig: Wig instance
        width: window half-width for local maximum detection (bp)
        height: minimum absolute signal value for summit calling.
                For mode='abs', this is the minimum |value|.
                For mode='pos', minimum value (positive peaks only).
                For mode='neg', minimum |value| for negative peaks only.
        mode: 'abs' = call peaks in both directions (|signal| >= height)
              'pos' = call only positive peaks (signal >= height)
              'neg' = call only negative peaks (signal <= -height)

    Returns:
        Summits instance
    """
    smts = Summits()
    step = wig.step
    half_w = max(1, width // (2 * step))
    back_w = half_w + 1

    for chrom in wig.data:
        arr = wig.data[chrom]
        n = arr.size

        if n < 2 * back_w:
            smts.data[chrom] = {'p': np.array([], dtype=int),
                                'v': np.array([], dtype=np.float64)}
            continue

        positions = []
        values = []

        i = half_w
        while i < n - back_w:
            window = arr[i - half_w:i + back_w]

            if mode == 'pos':
                v = window.max()
                if v < height:
                    i += back_w
                    continue
                if v == arr[i]:
                    # Find extent of flat top
                    start, end = i, i
                    while start > 0 and arr[start - 1] == v:
                        start -= 1
                    while end < n - 1 and arr[end + 1] == v:
                        end += 1
                    center = (start + end) // 2
                    positions.append(center * step)
                    values.append(v)
                    i = end + back_w
                else:
                    i += 1

            elif mode == 'neg':
                v = window.min()
                if v > -height:
                    i += back_w
                    continue
                if v == arr[i]:
                    start, end = i, i
                    while start > 0 and arr[start - 1] == v:
                        start -= 1
                    while end < n - 1 and arr[end + 1] == v:
                        end += 1
                    center = (start + end) // 2
                    positions.append(center * step)
                    values.append(v)
                    i = end + back_w
                else:
                    i += 1

            else:  # mode == 'abs'
                abs_window = np.abs(window)
                v_abs = abs_window.max()
                if v_abs < height:
                    i += back_w
                    continue
                if v_abs == abs(arr[i]):
                    start, end = i, i
                    target = arr[i]
                    while start > 0 and arr[start - 1] == target:
                        start -= 1
                    while end < n - 1 and arr[end + 1] == target:
                        end += 1
                    center = (start + end) // 2
                    positions.append(center * step)
                    values.append(arr[center])
                    i = end + back_w
                else:
                    i += 1

        smts.data[chrom] = {
            'p': np.array(positions, dtype=int),
            'v': np.array(values, dtype=np.float64)
        }
        print(f"  {chrom}: {len(positions)} summits called")

    return smts


# ---------------------------------------------------------------------------
# Edge detection: find position boundaries
# ---------------------------------------------------------------------------

def find_edges(wig, summits, distance=100, height=0.1):
    """
    Find start/end boundaries for each summit position.

    Searches outward from each summit until the signal crosses
    the height threshold, reaches a local minimum, or hits a
    neighboring summit boundary.

    Returns:
        List of dicts: [{chrom, start, end, summit_pos, summit_value, zscore, ...}, ...]
    """
    step = wig.step
    half_dist = distance // (2 * step)
    results = []

    for chrom in sorted(summits.data.keys()):
        if chrom not in wig.data:
            continue

        arr = wig.data[chrom]
        max_idx = arr.size - 1
        ps = summits.data[chrom]['p']
        vs = summits.data[chrom]['v']

        has_z = 'zscore' in summits.data[chrom]
        has_fuz = 's' in summits.data[chrom]

        for i in range(ps.size):
            pos_idx = ps[i] // step
            val = vs[i]

            # Determine threshold for this summit
            if val >= 0:
                thresh = height
            else:
                thresh = -height

            # Search left for start
            start_idx = int(pos_idx) - 1
            while start_idx > 0:
                if abs(start_idx - pos_idx) > half_dist:
                    break
                if i > 0 and start_idx * step <= ps[i - 1]:
                    start_idx += 1
                    break
                if val >= 0 and arr[start_idx] < thresh:
                    break
                if val < 0 and arr[start_idx] > thresh:
                    break
                # Local minimum check
                if abs(start_idx - pos_idx) > half_dist // 2:
                    lo = max(0, start_idx - 2)
                    hi = min(max_idx, start_idx + 3)
                    if val >= 0:
                        if arr[start_idx] == arr[lo:hi].min():
                            break
                    else:
                        if arr[start_idx] == arr[lo:hi].max():
                            break
                start_idx -= 1
            start_idx = max(0, start_idx)

            # Search right for end
            end_idx = int(pos_idx) + 1
            while end_idx < max_idx:
                if abs(end_idx - pos_idx) > half_dist:
                    break
                if i < ps.size - 1 and end_idx * step >= ps[i + 1]:
                    end_idx -= 1
                    break
                if val >= 0 and arr[end_idx] < thresh:
                    break
                if val < 0 and arr[end_idx] > thresh:
                    break
                if abs(end_idx - pos_idx) > half_dist // 2:
                    lo = max(0, end_idx - 2)
                    hi = min(max_idx + 1, end_idx + 3)
                    if val >= 0:
                        if arr[end_idx] == arr[lo:hi].min():
                            break
                    else:
                        if arr[end_idx] == arr[lo:hi].max():
                            break
                end_idx += 1
            end_idx = min(max_idx, end_idx)

            rec = {
                'chrom': chrom,
                'start': start_idx * step + 1,  # 1-based
                'end': end_idx * step + 1,
                'summit_pos': ps[i] + 1,  # 1-based
                'summit_value': val,
            }

            if has_z:
                rec['zscore'] = summits.data[chrom]['zscore'][i]
            if has_fuz:
                rec['fuzziness'] = summits.data[chrom]['s'][i]
                rec['fuz_log10p'] = summits.data[chrom]['ppos'][i]

            results.append(rec)

    return results


# ---------------------------------------------------------------------------
# Region calling: broad signal domains
# ---------------------------------------------------------------------------

def call_regions(wig, height=0.1, min_width=40, distance=100, mode='abs'):
    """
    Call broad signal regions (domains above/below threshold).

    Parameters:
        wig: Wig instance
        height: signal threshold (absolute value for mode='abs')
        min_width: minimum region width in bp
        distance: merge distance for neighboring regions
        mode: 'abs', 'pos', or 'neg'

    Returns:
        dict: regions[chrom] = [(start, end, summit_pos, summit_val, auc), ...]
    """
    step = wig.step
    min_bins = max(1, min_width // step)
    merge_bins = distance // step
    regions = {}

    for chrom in sorted(wig.data.keys()):
        arr = wig.data[chrom]

        if mode == 'pos':
            mask = arr >= height
        elif mode == 'neg':
            mask = arr <= -height
        else:  # abs
            mask = np.abs(arr) >= height

        # Find contiguous runs above threshold
        changes = np.diff(mask.astype(int))
        starts = np.where(changes == 1)[0] + 1
        ends = np.where(changes == -1)[0] + 1

        # Handle edge cases
        if mask[0]:
            starts = np.insert(starts, 0, 0)
        if mask[-1]:
            ends = np.append(ends, arr.size)

        # Filter by minimum width
        valid = (ends - starts) >= min_bins
        starts = starts[valid]
        ends = ends[valid]

        # Merge nearby regions
        if len(starts) > 1:
            merged_starts = [starts[0]]
            merged_ends = [ends[0]]
            for j in range(1, len(starts)):
                if starts[j] - merged_ends[-1] <= merge_bins:
                    merged_ends[-1] = ends[j]
                else:
                    merged_starts.append(starts[j])
                    merged_ends.append(ends[j])
            starts = np.array(merged_starts)
            ends = np.array(merged_ends)

        # Compute summit and AUC for each region
        chrom_regions = []
        for j in range(len(starts)):
            s, e = starts[j], ends[j]
            segment = arr[s:e]

            if mode == 'neg':
                smt_local = np.argmin(segment)
                smt_val = segment[smt_local]
            elif mode == 'pos':
                smt_local = np.argmax(segment)
                smt_val = segment[smt_local]
            else:
                abs_seg = np.abs(segment)
                smt_local = np.argmax(abs_seg)
                smt_val = segment[smt_local]

            smt_pos = (s + smt_local) * step
            auc = float(np.sum(segment)) * step

            chrom_regions.append((
                s * step,       # start bp
                e * step,       # end bp
                smt_val,        # summit value
                smt_pos,        # summit position
                auc             # area under curve
            ))

        regions[chrom] = chrom_regions
        print(f"  {chrom}: {len(chrom_regions)} regions")

    return regions


# ---------------------------------------------------------------------------
# Full position calling pipeline
# ---------------------------------------------------------------------------

def call_positions(wig, width=40, distance=100, height=0.1, mode='abs',
                   edge=True, fill_gap=False, fill_height=0.05,
                   compute_fuzziness=False, fuzziness_rd=None):
    """
    Full position calling pipeline: summit → merge → edges → output.

    Parameters:
        wig: Wig instance with signal data
        width: sliding window width for summit detection (bp)
        distance: minimum center-to-center distance between positions (bp)
        height: minimum signal threshold (absolute value for mode='abs')
        mode: 'abs' (both directions), 'pos' (positive only), 'neg' (negative only)
        edge: whether to search for position edges
        fill_gap: fill gaps between positions
        fill_height: minimum signal for gap filling
        compute_fuzziness: calculate fuzziness scores
        fuzziness_rd: half-region size for fuzziness (default: distance)

    Returns:
        (summits, positions): Summits instance and list of position dicts
    """
    if fuzziness_rd is None:
        fuzziness_rd = distance

    print("Calling summits...")
    smts = call_summits(wig, width=width, height=height, mode=mode)

    print("Merging nearby summits...")
    smts.merge(wig, distance=distance)

    if fill_gap:
        print("Filling gaps...")
        smts.fillgap(wig, height=fill_height, distance=distance)

    print("Computing z-scores...")
    smts.compute_zscores(wig)

    if compute_fuzziness:
        print("Computing fuzziness scores...")
        smts.compute_fuzziness(wig, rd=fuzziness_rd)

    if edge:
        print("Detecting position edges...")
        positions = find_edges(wig, smts, distance=distance, height=height)
    else:
        # No edge detection: use fixed flanking size
        flank = 74  # ~half nucleosome
        positions = []
        for chrom in sorted(smts.data.keys()):
            ps = smts.data[chrom]['p']
            vs = smts.data[chrom]['v']
            has_z = 'zscore' in smts.data[chrom]
            has_fuz = 's' in smts.data[chrom]
            for i in range(ps.size):
                rec = {
                    'chrom': chrom,
                    'start': max(1, ps[i] - flank + 1),
                    'end': ps[i] + flank + 1,
                    'summit_value': vs[i],
                    'summit_pos': ps[i] + 1,
                }
                if has_z:
                    rec['zscore'] = smts.data[chrom]['zscore'][i]
                if has_fuz:
                    rec['fuzziness'] = smts.data[chrom]['s'][i]
                    rec['fuz_log10p'] = smts.data[chrom]['ppos'][i]
                positions.append(rec)

    total = sum(smts.data[c]['p'].size for c in smts.data)
    print(f"Total positions called: {total}")

    return smts, positions


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def write_positions(positions, outfile):
    """Write position results to tab-delimited file."""
    if not positions:
        print("No positions to write.")
        return

    # Determine columns from first record
    cols = ['chrom', 'start', 'end', 'summit_value', 'summit_pos']
    extra = [k for k in positions[0] if k not in cols]
    cols.extend(sorted(extra))

    with open(outfile, 'w') as f:
        f.write('#' + '\t'.join(cols) + '\n')
        for rec in positions:
            vals = []
            for c in cols:
                v = rec.get(c, 'NA')
                if isinstance(v, float):
                    vals.append(f"{v:.6f}")
                else:
                    vals.append(str(v))
            f.write('\t'.join(vals) + '\n')

    print(f"Wrote {len(positions)} positions to {outfile}")


def write_regions(regions, outfile):
    """Write region results to tab-delimited file."""
    with open(outfile, 'w') as f:
        f.write('#chrom\tstart\tend\tsummit_value\tsummit_pos\tauc\twidth\n')
        total = 0
        for chrom in sorted(regions.keys()):
            for (start, end, smt_pos, smt_val, auc) in regions[chrom]:
                width = end - start
                f.write(f"{chrom}\t{start}\t{end}\t{smt_val:.6f}\t{smt_pos}\t"
                        f"{auc:.2f}\t{width}\n")
                total += 1
    print(f"Wrote {total} regions to {outfile}")


# ---------------------------------------------------------------------------
# Differential analysis between two wig files
# ---------------------------------------------------------------------------

def differential_positions(wig1, wig2, positions):
    """
    Compute differential signal between two Wig instances at each position.

    Adds 'value_1', 'value_2', 'diff', and 'diff_zscore' to each position record.
    Uses a Welch's t-test on the signal within each position's boundaries.
    """
    for rec in positions:
        chrom = rec['chrom']
        s = (rec['start'] - 1) // wig1.step
        e = rec['end'] // wig1.step

        if chrom in wig1.data and chrom in wig2.data:
            s1 = min(s, wig1.data[chrom].size)
            e1 = min(e, wig1.data[chrom].size)
            s2 = min(s, wig2.data[chrom].size)
            e2 = min(e, wig2.data[chrom].size)

            region1 = wig1.data[chrom][s1:e1]
            region2 = wig2.data[chrom][s2:e2]

            rec['value_1'] = float(np.mean(region1)) if region1.size > 0 else 0.0
            rec['value_2'] = float(np.mean(region2)) if region2.size > 0 else 0.0
            rec['diff'] = rec['value_1'] - rec['value_2']

            if region1.size > 1 and region2.size > 1:
                try:
                    t_stat, p_val = stats.ttest_ind(region1, region2, equal_var=False)
                    if p_val > 0:
                        rec['diff_log10p'] = -np.log10(p_val)
                    else:
                        rec['diff_log10p'] = 300.0
                except:
                    rec['diff_log10p'] = 0.0
            else:
                rec['diff_log10p'] = 0.0
        else:
            rec['value_1'] = 0.0
            rec['value_2'] = 0.0
            rec['diff'] = 0.0
            rec['diff_log10p'] = 0.0

    return positions


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='signal_peaks: Position/peak/region calling for continuous '
                    'normalized coverage signals [-1, 1].',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('command', choices=['dpos', 'dregion', 'dtriple'],
                        help="'dpos': call positions, 'dregion': call regions, "
                             "'dtriple': call both")
    parser.add_argument('wig_file', help='Path to input bigWig (.bw) or wiggle (.wig) file')

    # Signal parameters
    parser.add_argument('-q', '--height', type=float, default=0.1,
                        help='Minimum absolute signal threshold for calling')
    parser.add_argument('--mode', choices=['abs', 'pos', 'neg'], default='abs',
                        help="'abs': bidirectional, 'pos': positive only, "
                             "'neg': negative only")
    parser.add_argument('-a', '--step', type=int, default=0,
                        help='Step/bin size in bp (0=auto: 1bp for bigWig, '
                             'detect from header for wig)')

    # Position calling parameters
    parser.add_argument('-jw', '--width', type=int, default=40,
                        help='Sliding window width for summit detection (bp)')
    parser.add_argument('-jd', '--distance', type=int, default=100,
                        help='Minimum center-to-center distance between positions')
    parser.add_argument('-e', '--edge', type=int, default=1,
                        help='Detect position edges (1=yes, 0=no)')
    parser.add_argument('-g', '--gapfill', type=int, default=0,
                        help='Fill gaps between positions (1=yes, 0=no)')
    parser.add_argument('--fuzziness', type=int, default=0,
                        help='Compute fuzziness/positioning scores (1=yes, 0=no)')

    # Region calling parameters
    parser.add_argument('-rw', '--region_width', type=int, default=40,
                        help='Minimum region width (bp)')
    parser.add_argument('-rd', '--region_dis', type=int, default=3000,
                        help='Merge distance for neighboring regions (bp)')

    # Smoothing
    parser.add_argument('-z', '--smooth_width', type=int, default=0,
                        help='Smoothing window width (bp), 0=no smoothing')

    # Differential analysis
    parser.add_argument('--wig2', default=None,
                        help='Second bigWig/wiggle file for differential analysis')

    # Output
    parser.add_argument('-o', '--out', default='result',
                        help='Output directory name')

    args = parser.parse_args()


if __name__ == '__main__':
    main()
