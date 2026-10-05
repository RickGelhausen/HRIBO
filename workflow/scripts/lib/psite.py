"""
Ribosomal-site offset estimation and read length scoring at translation boundaries.

A translation initiation site caller such as ORFBounder needs to know which read
lengths carry a usable initiation signal, and by how much to shift each of them
so that a read is attributed to the codon its ribosome actually occupied. This
module derives both from the metagene profiles the workflow already computes.

The offset is the distance from the mapped read end to the calibrated ribosomal
site: the P-site at initiation and the A-site at termination. All coordinates
run in the direction of translation, with zero at the first nucleotide of the
selected start or stop codon. For a stop-codon A-site calibration, the P-site is
three nucleotides upstream; frame and periodicity evidence is taken from the
coding body before the stop, excluding the terminal peak.

A deliberate caveat runs through the scoring: three-nucleotide periodicity is
frequently weak or absent in bacterial Ribo-seq, far more so than in eukaryotes.
Periodicity is therefore reported, but a read length is not rejected
for lacking it, and the confidence attached to a recommendation says which of
the two lines of evidence it rests on.

Author: Rick Gelhausen
"""

from __future__ import annotations

import functools
from dataclasses import dataclass, field

import numpy as np

@dataclass(frozen=True)
class ReadEnd:
    """The geometry of measuring offsets from one end of a read.

    Which end carries the cleaner signal is organism and protocol dependent: in
    many bacteria nuclease digestion defines the 3' end more sharply than the
    5' end, while the reverse holds elsewhere. Both are therefore first class
    here rather than one being assumed.

    sign             +1 when the mapped end lies upstream of the site (5'),
                     -1 when it lies downstream (3')
    search_window    where the boundary peak is looked for, relative to the
                     selected codon. A 5' end peaks upstream of it, a 3' end
                     downstream, roughly a read length away
    plausible        offsets outside this range indicate noise rather than a
                     genuine boundary peak
    """

    name: str
    sign: int
    search_window: tuple[int, int]
    plausible: tuple[int, int]

    def offset_from_peak(self, peak_position: int) -> int:
        """Convert the peak's position into a positive distance to the site."""
        return -self.sign * peak_position

    def site_shift(self, offset: int) -> int:
        """Move a profile onto the site used to calibrate its offset."""
        return self.sign * offset

    def psite_shift(self, offset: int) -> int:
        """Historical alias for a start-codon P-site calibration's shift."""
        return self.site_shift(offset)

    def is_plausible(self, offset: int) -> bool:
        return self.plausible[0] <= offset <= self.plausible[1]


READ_ENDS = {
    # A ribosome with its P-site on the start codon puts the read's 5' end about
    # 12 nt upstream of it.
    "fiveprime": ReadEnd("fiveprime", 1, (-25, 5), (5, 20)),
    # The same ribosome puts the read's 3' end downstream of the start codon, by
    # the read length minus the 5' offset, so the plausible range is wider.
    "threeprime": ReadEnd("threeprime", -1, (-5, 32), (5, 28)),
}

DEFAULT_READ_END = READ_ENDS["fiveprime"]

# The A-site lies three nucleotides downstream of the P-site. The corresponding
# distances from a 5' end grow by three; distances from a 3' end shrink by three.
# Search windows retain a margin to diagnose an implausible maximum, rather
# than choosing a smaller incidental peak inside the accepted offset range.
STOP_READ_ENDS = {
    "fiveprime": ReadEnd("fiveprime", 1, (-28, 5), (8, 23)),
    "threeprime": ReadEnd("threeprime", -1, (-5, 29), (2, 25)),
}


def _validate_anchor(anchor: str) -> None:
    if anchor not in {"start", "stop"}:
        raise ValueError("anchor must be 'start' or 'stop'")


def _read_end_for_anchor(
    read_end: ReadEnd, anchor: str, read_length: int | None = None
) -> ReadEnd:
    _validate_anchor(anchor)
    if anchor == "start":
        return read_end
    if read_length is not None and read_length >= 50 and read_end.name == "fiveprime":
        # Long termination footprints can span a queued ribosome pair. The
        # terminating ribosome remains near the 3' end; mirror its 3' A-site
        # distances to search for that same site from the distant 5' end.
        low = read_length - 1 - STOP_READ_ENDS["threeprime"].plausible[1]
        high = read_length - 1 - STOP_READ_ENDS["threeprime"].plausible[0]
        return ReadEnd("fiveprime", 1, (-high - 5, -low + 5), (low, high))
    # Preserve explicitly supplied geometry, including the length-dependent
    # disome window chosen by the caller above.
    return STOP_READ_ENDS[read_end.name] if read_end == READ_ENDS[read_end.name] else read_end


def read_end_for_anchor(
    name: str, anchor: str = "start", read_length: int | None = None
) -> ReadEnd:
    """Return read-end geometry for a start/P-site or stop/A-site calibration.

    Stop footprints of at least 50 nt use a 5' window consistent with a queued
    ribosome pair and the same terminating-ribosome 3' offset bounds. This does
    not establish whether a long footprint actually represents a disome.
    """
    return _read_end_for_anchor(READ_ENDS[name], anchor, read_length)


def _site_for_anchor(anchor: str) -> str:
    _validate_anchor(anchor)
    return "A" if anchor == "stop" else "P"


def _signal_for_anchor(anchor: str) -> str:
    _validate_anchor(anchor)
    return "termination" if anchor == "stop" else "initiation"


def psite_offset(
    offset: int, read_end: ReadEnd = DEFAULT_READ_END, anchor: str = "start"
) -> int:
    """Convert a calibrated offset into the corresponding P-site distance.

    Stop profiles calibrate the A-site: a 5' P-site distance is A - 3, and a 3'
    P-site distance is A + 3. Start profiles already calibrate the P-site.
    """
    _validate_anchor(anchor)
    return offset - read_end.sign * 3 if anchor == "stop" else offset


# A read length carrying less than this share of the library is too sparse to
# base an offset on, however clean its profile looks.
MIN_ABUNDANCE_FRACTION = 0.005
MIN_BACKGROUND_POSITIONS = 10
MIN_FRAME_READS = 30

# A read length joins the recommended set only if it improves the pooled peak by
# at least this much, relatively. Accepting any improvement at all lets a read
# length with no real signal in on a rounding difference.
MIN_RELATIVE_GAIN = 0.02

# Two read ends whose pooled peaks are within this factor of each other are
# treated as equally sharp. This is the normal case rather than the exception:
# for a read of fixed length the 5' and 3' positions are rigidly linked, so the
# two per-read-length profiles are shifted copies of one another and are equally
# sharp by construction. Sharpness only separates the ends when one of them
# fails outright.
SHARPNESS_TIE_FACTOR = 1.25

# Frame bias differences smaller than this are not meaningful.
FRAME_TIE_MARGIN = 0.05

# Region treated as background when judging how far the peak stands out. Taken
# well upstream of the start codon, where little initiation signal is expected.
BACKGROUND_WINDOW = (-100, -30)
STOP_BACKGROUND_WINDOW = (-100, -40)


def _background_window(anchor: str) -> tuple[int, int]:
    return STOP_BACKGROUND_WINDOW if anchor == "stop" else BACKGROUND_WINDOW


def _body_mask(psite_positions: np.ndarray, body_length: int, anchor: str) -> np.ndarray:
    """Exclude the boundary peak and select coding-body positions."""
    if anchor == "stop":
        return (psite_positions >= -15 - body_length) & (psite_positions < -15)
    return (psite_positions >= 15) & (psite_positions < 15 + body_length)


def _observed_site_mask(coordinates, site_positions, valid_mask=None):
    """Keep observations that survive alignment onto the configured output axis."""
    coordinates = np.asarray(coordinates)
    if not coordinates.size:
        return np.zeros(0, dtype=bool)
    valid = (site_positions >= coordinates[0]) & (site_positions <= coordinates[-1])
    if valid_mask is not None:
        valid &= np.asarray(valid_mask, dtype=bool)
    return valid


def background_statistics(profile, site_coordinates, anchor="start", valid_mask=None):
    """Measure the same observed site-coordinate background for every metric.

    Zero-filled bins introduced by aligning profiles are not observations.
    A zero median falls back to the background mean, then the observed profile
    mean for a genuinely empty baseline, keeping enrichment finite.
    """
    _validate_anchor(anchor)
    profile = np.asarray(profile, dtype=float)
    site_coordinates = np.asarray(site_coordinates)
    valid = np.ones(profile.size, dtype=bool) if valid_mask is None else np.asarray(valid_mask, dtype=bool)
    values = profile[valid & _window_slice(site_coordinates, *_background_window(anchor))]
    if values.size < MIN_BACKGROUND_POSITIONS:
        return 0.0, 0.0, 0.0, 0.0, int(values.size)
    median = float(np.median(values))
    mean = float(values.mean())
    sd = float(values.std())
    reference = median if median > 0 else mean
    if reference <= 0:
        reference = float(profile[valid].mean()) if valid.any() else 0.0
    return median, mean, sd, reference, int(values.size)


# How far above the background a peak must stand, in background standard
# deviations. A ratio alone is not enough: in a shallow profile the largest
# value of pure Poisson noise routinely reaches three or four times the median,
# so a ratio test accepts noise as an initiation peak. A real initiation peak
# clears this threshold by a wide margin.
MIN_PEAK_Z = 5.0

# The peak must also be this many times the background, which keeps very deep
# but flat profiles from passing on the z-score alone.
MIN_PEAK_RATIO = 2.0


@dataclass
class ReadLengthScore:
    """Everything known about one read length, and what it is worth."""

    read_length: int
    total_reads: int
    abundance: float

    offset: int | None
    peak_height: float
    background: float
    sharpness: float
    z_score: float

    frame_fractions: tuple[float, float, float]
    frame_bias: float
    periodicity: float

    usable: bool
    reasons: list[str] = field(default_factory=list)
    anchor: str = "start"
    site: str = "P"
    background_reference: float = 0.0
    background_positions: int = 0
    frame_reads: float | None = None

    @property
    def score(self) -> float:
        """Peak/background enrichment, without a ceiling that erases differences."""
        return float(self.sharpness) if self.usable else 0.0


def _window_slice(coordinates: np.ndarray, low: int, high: int) -> np.ndarray:
    return (coordinates >= low) & (coordinates <= high)


@dataclass
class PeakEstimate:
    """The boundary peak of one profile, and how convincing it is."""

    offset: int | None
    peak_height: float
    background: float
    background_sd: float
    sharpness: float
    z_score: float
    background_reference: float = 0.0
    background_positions: int = 0

    @property
    def is_significant(self) -> bool:
        return (
            self.offset is not None
            and self.z_score >= MIN_PEAK_Z
            and self.sharpness >= MIN_PEAK_RATIO
        )


def estimate_offset(
    profile: np.ndarray,
    coordinates: np.ndarray,
    read_end: ReadEnd = DEFAULT_READ_END,
    search_window: tuple[int, int] | None = None,
    anchor: str = "start",
) -> PeakEstimate:
    """Locate a start/P-site or stop/A-site peak and assess its evidence.

    Coordinates must be relative to the first base of the chosen codon. The
    offset is the positive distance from the mapped read end to that codon, and
    is None when the search-window maximum falls at an implausible offset.
    """
    read_end = _read_end_for_anchor(read_end, anchor)
    profile = np.asarray(profile, dtype=float)
    coordinates = np.asarray(coordinates)
    search_window = search_window or read_end.search_window

    search = _window_slice(coordinates, *search_window)
    if not search.any() or profile[search].sum() == 0:
        return PeakEstimate(None, 0.0, 0.0, 0.0, 0.0, 0.0)

    search_positions = coordinates[search]
    search_values = profile[search]
    peak_position = int(search_positions[int(np.argmax(search_values))])
    peak_height = float(search_values.max())

    offset = read_end.offset_from_peak(peak_position)
    site_coordinates = coordinates + read_end.site_shift(offset)
    background, background_mean, background_sd, reference, background_positions = background_statistics(
        profile, site_coordinates, anchor, _observed_site_mask(coordinates, site_coordinates)
    )
    sharpness = float(peak_height / reference) if reference > 0 else 0.0

    # Counts are Poisson-like, so the variance is at least the mean; use that as
    # a floor to stop an unusually smooth background inflating the z-score.
    spread = max(background_sd, np.sqrt(max(background_mean, 0.0)), 1.0)
    z_score = float((peak_height - background_mean) / spread)

    if not read_end.is_plausible(offset) or background_positions < MIN_BACKGROUND_POSITIONS:
        offset = None

    return PeakEstimate(offset, peak_height, background, background_sd, sharpness, z_score,
                        reference, background_positions)


def frame_fractions(
    profile: np.ndarray,
    coordinates: np.ndarray,
    offset: int,
    read_end: ReadEnd = DEFAULT_READ_END,
    body_length: int = 90,
    anchor: str = "start",
    valid_mask: np.ndarray | None = None,
) -> tuple[float, float, float]:
    """Share of P-sites falling in each reading frame inside the ORF body.

    Frame 0 is the frame of the first base of the selected boundary codon. The
    body lies downstream of a start or upstream of a stop, excluding the peak.
    Stop-calibrated A-site offsets are converted to P-site positions first.
    """
    _validate_anchor(anchor)
    profile = np.asarray(profile, dtype=float)
    coordinates = np.asarray(coordinates)

    psite_positions = coordinates + read_end.site_shift(psite_offset(offset, read_end, anchor))
    site_positions = coordinates + read_end.site_shift(offset)
    body = _body_mask(psite_positions, body_length, anchor) & _observed_site_mask(
        coordinates, site_positions, valid_mask
    )
    if not body.any():
        return (0.0, 0.0, 0.0)

    counts = np.zeros(3)
    for frame in range(3):
        selected = body & (psite_positions % 3 == frame)
        counts[frame] = profile[selected].sum()

    total = counts.sum()
    if total <= 0:
        return (0.0, 0.0, 0.0)
    return tuple(float(c / total) for c in counts)


def frame_read_count(profile, coordinates, offset, read_end=DEFAULT_READ_END,
                     body_length=90, anchor="start", valid_mask=None):
    """Observed coding-body read ends underlying the frame fractions."""
    _validate_anchor(anchor)
    coordinates = np.asarray(coordinates)
    positions = coordinates + read_end.site_shift(psite_offset(offset, read_end, anchor))
    site_positions = coordinates + read_end.site_shift(offset)
    body = _body_mask(positions, body_length, anchor) & _observed_site_mask(
        coordinates, site_positions, valid_mask
    )
    return float(np.asarray(profile)[body].sum())


def _supported_frame_zero(fractions, frame_reads=None):
    return fractions[0] if frame_reads is None or frame_reads >= MIN_FRAME_READS else 0.0


def periodicity_score(
    profile: np.ndarray,
    coordinates: np.ndarray,
    offset: int,
    read_end: ReadEnd = DEFAULT_READ_END,
    body_length: int = 90,
    anchor: str = "start",
    valid_mask: np.ndarray | None = None,
) -> float:
    """Strength of the 3-nt component of the elongation signal, in [0, 1].

    Computed as the share of the spectrum's power that sits at a period of three
    nucleotides, which is less sensitive to a single dominant codon than the
    frame fractions are.
    """
    _validate_anchor(anchor)
    profile = np.asarray(profile, dtype=float)
    coordinates = np.asarray(coordinates)

    psite_positions = coordinates + read_end.site_shift(psite_offset(offset, read_end, anchor))
    site_positions = coordinates + read_end.site_shift(offset)
    body = _body_mask(psite_positions, body_length, anchor) & _observed_site_mask(
        coordinates, site_positions, valid_mask
    )
    values = profile[body]
    if values.size < 9 or values.sum() <= 0:
        return 0.0

    values = values - values.mean()
    if not np.any(values):
        return 0.0

    spectrum = np.abs(np.fft.rfft(values)) ** 2
    frequencies = np.fft.rfftfreq(values.size)
    if spectrum[1:].sum() <= 0:
        return 0.0

    # Bin closest to one cycle per three nucleotides.
    target = int(np.argmin(np.abs(frequencies - 1 / 3)))
    if target == 0:
        return 0.0
    return float(spectrum[target] / spectrum[1:].sum())


def score_read_lengths(
    start_profiles: dict[int, np.ndarray],
    coordinates: np.ndarray,
    read_totals: dict[int, int] | None = None,
    read_end: ReadEnd = DEFAULT_READ_END,
    anchor: str = "start",
) -> list[ReadLengthScore]:
    """Score every read length in the selected boundary's metagene profile.

    The historical ``start_profiles`` argument also accepts stop profiles when
    ``anchor='stop'``; coordinates must then use the first stop base as zero.
    """
    read_end = _read_end_for_anchor(read_end, anchor)
    signal = _signal_for_anchor(anchor)
    if read_totals is None:
        read_totals = {
            length: int(np.asarray(profile).sum()) for length, profile in start_profiles.items()
        }
    library_total = sum(read_totals.values()) or 1

    scores: list[ReadLengthScore] = []
    for read_length in sorted(start_profiles):
        profile = np.asarray(start_profiles[read_length], dtype=float)
        total = int(read_totals.get(read_length, 0))
        abundance = total / library_total

        length_read_end = _read_end_for_anchor(read_end, anchor, read_length)
        estimate = estimate_offset(profile, coordinates, length_read_end, anchor=anchor)
        offset = estimate.offset

        reasons: list[str] = []
        usable = True

        if abundance < MIN_ABUNDANCE_FRACTION:
            usable = False
            reasons.append(
                f"carries only {abundance:.2%} of the library, below the {MIN_ABUNDANCE_FRACTION:.1%} minimum"
            )

        if offset is not None and offset >= read_length:
            usable = False
            reasons.append(
                f"the estimated {_site_for_anchor(anchor)}-site offset {offset} nt lies "
                f"outside a {read_length}-nt read (offset must be less than read length)"
            )
            offset = None

        if estimate.offset is None:
            usable = False
            if estimate.peak_height > 0 and estimate.background_positions < MIN_BACKGROUND_POSITIONS:
                reasons.append(
                    f"insufficient observed background ({estimate.background_positions} positions; "
                    f"at least {MIN_BACKGROUND_POSITIONS} are required); widen the profile window"
                )
            else:
                reasons.append(
                    f"no {signal} peak within "
                    f"{length_read_end.plausible[0]}-{length_read_end.plausible[1]} nt of the {anchor} codon"
                )

        if offset is None:
            fractions = (0.0, 0.0, 0.0)
            periodicity = 0.0
            frame_reads = 0.0
        else:
            fractions = frame_fractions(profile, coordinates, offset, read_end, anchor=anchor)
            periodicity = periodicity_score(profile, coordinates, offset, read_end, anchor=anchor)
            frame_reads = frame_read_count(profile, coordinates, offset, read_end, anchor=anchor)
            if not estimate.is_significant:
                usable = False
                reasons.append(
                    f"the peak at offset {estimate.offset} is not distinguishable from noise "
                    f"({estimate.sharpness:.1f}x background, {estimate.z_score:.1f} standard "
                    f"deviations; at least {MIN_PEAK_RATIO:.0f}x and {MIN_PEAK_Z:.0f} are required)"
                )

        scores.append(
            ReadLengthScore(
                read_length=read_length,
                total_reads=total,
                abundance=abundance,
                offset=offset,
                peak_height=estimate.peak_height,
                background=estimate.background,
                sharpness=estimate.sharpness,
                z_score=estimate.z_score,
                frame_fractions=fractions,
                frame_bias=max(fractions) if any(fractions) else 0.0,
                periodicity=periodicity,
                usable=usable,
                reasons=reasons,
                anchor=anchor,
                site=_site_for_anchor(anchor),
                background_reference=estimate.background_reference,
                background_positions=estimate.background_positions,
                frame_reads=frame_reads,
            )
        )

    return scores


# --------------------------------------------------------------------------
# Combining read lengths
# --------------------------------------------------------------------------


@dataclass
class Recommendation:
    """Suggested boundary profiling parameters, and how much to trust them."""

    read_lengths: list[int]
    offsets: dict[int, int]
    read_end: str

    sharpness: float
    frame_fractions: tuple[float, float, float]
    frame_bias: float
    periodicity: float
    covered_fraction: float

    confidence: str
    rationale: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    anchor: str = "start"
    site: str = "P"
    frame_reads: float | None = None

    @property
    def has_recommendation(self) -> bool:
        return bool(self.read_lengths)

    def offset_table(self) -> list[tuple[int, int]]:
        return [(length, self.offsets[length]) for length in sorted(self.read_lengths)]


def shift_profile(profile: np.ndarray, offset: int) -> np.ndarray:
    """Shift a profile from read-end coordinates onto a calibrated site.

    A read end at c is attributed to c + offset. The profile moves ``offset``
    positions to the right, zero filled. The offset is signed here.
    """
    profile = np.asarray(profile, dtype=float)
    if offset == 0:
        return profile.copy()
    shifted = np.zeros_like(profile)
    if abs(offset) >= profile.size:
        return shifted
    if offset > 0:
        shifted[offset:] = profile[: profile.size - offset]
    else:
        shifted[: profile.size + offset] = profile[-offset:]
    return shifted


def pool_profiles(
    profiles: dict[int, np.ndarray],
    offsets: dict[int, int],
    read_lengths: list[int],
    read_end: ReadEnd = DEFAULT_READ_END,
) -> np.ndarray:
    """Sum the offset-corrected profiles of several read lengths."""
    pooled = None
    for read_length in read_lengths:
        shifted = shift_profile(
            profiles[read_length], read_end.site_shift(offsets[read_length])
        )
        pooled = shifted if pooled is None else pooled + shifted
    if pooled is None:
        raise ValueError("pool_profiles requires at least one read length")
    return pooled


def _evaluate_combination(
    profiles: dict[int, np.ndarray],
    coordinates: np.ndarray,
    offsets: dict[int, int],
    read_lengths: list[int],
    read_end: ReadEnd = DEFAULT_READ_END,
    anchor: str = "start",
) -> tuple[float, tuple[float, float, float], float]:
    """Sharpness, frame fractions and periodicity of a pooled set of lengths."""
    pooled = pool_profiles(profiles, offsets, read_lengths, read_end)

    # Offset correction aligns the selected boundary's P- or A-site at zero.
    at_boundary = pooled[coordinates == 0]
    peak = float(at_boundary[0]) if at_boundary.size else 0.0

    valid = _combination_valid_mask(pooled.size, offsets, read_lengths, read_end)
    _, _, _, reference, _ = background_statistics(pooled, coordinates, anchor, valid)
    sharpness = float(peak / reference) if reference > 0 else 0.0

    # A-site-aligned termination profiles still need the 3-nt P-site correction
    # when evaluating coding-body evidence.
    fractions = frame_fractions(pooled, coordinates, 0, read_end, anchor=anchor, valid_mask=valid)
    periodicity = periodicity_score(pooled, coordinates, 0, read_end, anchor=anchor, valid_mask=valid)
    return sharpness, fractions, periodicity


def _combination_valid_mask(size, offsets, read_lengths, read_end):
    """Positions observed for every member of an aligned pool."""
    valid = np.ones(size, dtype=bool)
    for length in read_lengths:
        valid &= shift_profile(np.ones(size), read_end.site_shift(offsets[length])).astype(bool)
    return valid


def _combination_frame_read_count(profiles, coordinates, offsets, read_lengths, read_end, anchor):
    pooled = pool_profiles(profiles, offsets, read_lengths, read_end)
    valid = _combination_valid_mask(pooled.size, offsets, read_lengths, read_end)
    return frame_read_count(pooled, coordinates, 0, read_end, anchor=anchor, valid_mask=valid)


def recommend_read_lengths(
    scores: list[ReadLengthScore],
    start_profiles: dict[int, np.ndarray],
    coordinates: np.ndarray,
    read_end: ReadEnd = DEFAULT_READ_END,
    anchor: str = "start",
) -> Recommendation:
    """Choose read lengths and calibrated-site offsets for boundary profiling.

    Start with the strongest single-length enrichment. At each step choose the
    remaining length with the largest gain in that same pooled metric. Frame 0
    support and read abundance only break enrichment ties. A marginal length
    must not dilute a clear peak merely because its coding-body frame is strong.
    """
    read_end = _read_end_for_anchor(read_end, anchor)
    signal = _signal_for_anchor(anchor)
    coordinates = np.asarray(coordinates)
    usable = [s for s in scores if s.usable]

    if not usable:
        return _no_recommendation(scores, read_end, anchor)

    offsets = {s.read_length: s.offset for s in usable}
    total_reads = sum(s.total_reads for s in scores) or 1

    single_metrics = {
        score.read_length: _evaluate_combination(
            start_profiles, coordinates, offsets, [score.read_length], read_end, anchor
        ) for score in usable
    }
    usable.sort(key=lambda s: (single_metrics[s.read_length][0],
                              _supported_frame_zero(s.frame_fractions, s.frame_reads),
                              s.total_reads, -s.read_length), reverse=True)
    selected = [usable[0].read_length]
    best_sharpness, best_fractions, best_periodicity = single_metrics[selected[0]]

    rationale = [
        f"read length {usable[0].read_length} has the strongest {signal} peak "
        f"({best_sharpness:.1f}x background at offset {usable[0].offset}); "
        "reading-frame bias does not override peak enrichment"
    ]

    remaining = usable[1:]
    while remaining:
        trials = [(candidate, _evaluate_combination(
            start_profiles, coordinates, offsets, selected + [candidate.read_length], read_end, anchor
        ), _combination_frame_read_count(
            start_profiles, coordinates, offsets, selected + [candidate.read_length], read_end, anchor
        )) for candidate in remaining]
        candidate, metrics, _ = max(trials, key=lambda trial: (
            trial[1][0], _supported_frame_zero(trial[1][1], trial[2]),
            trial[0].total_reads, -trial[0].read_length
        ))
        sharpness, fractions, periodicity = metrics
        if sharpness < best_sharpness * (1 + MIN_RELATIVE_GAIN):
            for excluded, (trial_sharpness, _, _), _ in trials:
                rationale.append(
                    f"read length {excluded.read_length} was left out: it moves the pooled "
                    f"peak from {best_sharpness:.1f}x to {trial_sharpness:.1f}x, short of the "
                    f"{MIN_RELATIVE_GAIN:.0%} gain required to include it"
                )
            break
        selected.append(candidate.read_length)
        remaining.remove(candidate)
        best_sharpness, best_fractions, best_periodicity = metrics
        rationale.append(
            f"adding read length {candidate.read_length} (offset {candidate.offset}) "
            f"raises the pooled peak to {sharpness:.1f}x"
        )

    covered = sum(s.total_reads for s in scores if s.read_length in selected) / total_reads
    frame_reads = _combination_frame_read_count(
        start_profiles, coordinates, offsets, selected, read_end, anchor
    )
    confidence, warnings = _assess_confidence(
        best_sharpness, best_fractions, best_periodicity, covered, anchor, frame_reads
    )

    return Recommendation(
        read_lengths=sorted(selected),
        offsets={length: offsets[length] for length in selected},
        read_end=read_end.name,
        sharpness=best_sharpness,
        frame_fractions=best_fractions,
        frame_bias=max(best_fractions) if any(best_fractions) else 0.0,
        periodicity=best_periodicity,
        covered_fraction=covered,
        confidence=confidence,
        rationale=rationale,
        warnings=warnings,
        anchor=anchor,
        site=_site_for_anchor(anchor),
        frame_reads=frame_reads,
    )


def _no_recommendation(
    scores: list[ReadLengthScore], read_end: ReadEnd = DEFAULT_READ_END,
    anchor: str = "start",
) -> Recommendation:
    """Explain why no read length is usable instead of inventing a setup."""
    signal = _signal_for_anchor(anchor)
    target = "TTS profiling" if anchor == "stop" else "TIS caller"
    warnings = [
        f"No read length shows a usable translation {signal} signal, so no {target} "
        "setup is suggested."
    ]
    if not scores:
        warnings.append("No reads were found in the metagene window at all.")
    else:
        best = max(scores, key=lambda s: s.sharpness, default=None)
        if best is not None and best.sharpness > 0 and (
            best.sharpness < MIN_PEAK_RATIO or best.z_score < MIN_PEAK_Z
        ):
            warnings.append(
                f"The closest candidate was read length {best.read_length} at "
                f"{best.sharpness:.1f}x background ({best.z_score:.1f} standard deviations), "
                f"below the required {MIN_PEAK_RATIO:.0f}x and {MIN_PEAK_Z:.0f}."
            )
        elif best is not None and best.reasons:
            warnings.append(
                f"The closest candidate was read length {best.read_length}: "
                + "; ".join(best.reasons)
                + "."
            )
        warnings.append(
            "Common causes: the library is RNA-seq rather than Ribo-seq, the reads are "
            f"dominated by rRNA or tRNA, or the annotated {anchor} codons are inaccurate."
        )
    if anchor == "stop":
        warnings.append(_termination_caveat())
    return Recommendation(
        read_lengths=[],
        offsets={},
        read_end=read_end.name,
        sharpness=0.0,
        frame_fractions=(0.0, 0.0, 0.0),
        frame_bias=0.0,
        periodicity=0.0,
        covered_fraction=0.0,
        confidence="none",
        rationale=[],
        warnings=warnings,
        anchor=anchor,
        site=_site_for_anchor(anchor),
    )


def _assess_confidence(
    sharpness: float,
    fractions: tuple[float, float, float],
    periodicity: float,
    covered: float,
    anchor: str = "start",
    frame_reads: float | None = None,
) -> tuple[str, list[str]]:
    """Grade a recommendation, and say which evidence it rests on."""
    frame_bias = max(fractions) if any(fractions) else 0.0
    in_frame = fractions[0] if any(fractions) else 0.0
    warnings: list[str] = []
    signal = _signal_for_anchor(anchor)

    if sharpness >= 5.0 and _supported_frame_zero(fractions, frame_reads) >= 0.5:
        confidence = "high"
    elif sharpness >= 5.0:
        confidence = "medium"
        if in_frame < 0.5:
            warnings.append(
                f"The {signal} peak is clear ({sharpness:.1f}x) but reading frame bias in frame 0 is "
                f"weak ({in_frame:.0%}). This is common in bacterial "
                f"Ribo-seq; the offsets rest on the {anchor}-codon peak alone."
            )
    elif sharpness >= 3.0:
        confidence = "medium"
    else:
        confidence = "low"
        warnings.append(
            f"The pooled {signal} peak is only {sharpness:.1f}x background. Treat the "
            "offsets as provisional and inspect the metagene plots before relying on them."
        )

    if frame_bias >= 0.5 and in_frame < frame_bias:
        dominant = int(np.argmax(fractions))
        warnings.append(
            f"The dominant reading frame is {dominant}, not 0. The estimated offsets may be "
            f"off by {dominant} nt, or the annotated {anchor} codons may be shifted."
        )

    if frame_reads is not None and frame_reads < MIN_FRAME_READS:
        warnings.append(
            f"Only {frame_reads:.0f} coding-body read ends support the frame fractions; "
            f"at least {MIN_FRAME_READS} are required to use frame bias for high confidence."
        )

    if covered < 0.3:
        warnings.append(
            f"The selected read lengths cover only {covered:.0%} of the library. Check the "
            "read length distribution for a broader usable range."
        )

    if periodicity < 0.15 and confidence != "low":
        warnings.append(
            f"Three-nucleotide periodicity is weak ({periodicity:.0%} of spectral power). "
            f"Sub-codon assignment will be unreliable even though the {signal} signal is good."
        )

    if anchor == "stop":
        warnings.append(_termination_caveat())

    return confidence, warnings


def _termination_caveat() -> str:
    return (
        "Termination peaks calibrate the stop codon in the A-site; the P-site is 3 nt "
        "upstream. Queued ribosomes, long footprints, and stop-codon readthrough can "
        "produce additional peaks; inspect the stop metagene before using the offsets."
    )


# --------------------------------------------------------------------------
# Choosing between the read ends
# --------------------------------------------------------------------------


@dataclass
class EndComparison:
    """The analysis for one read end, kept so both can be reported."""

    read_end: str
    scores: list[ReadLengthScore]
    recommendation: Recommendation

    @property
    def offset_spread(self) -> int | None:
        """How much the estimated offset varies across usable read lengths.

        This, rather than sharpness, is what actually distinguishes the two ends.
        A protocol that defines one end precisely produces the same offset at
        every read length from that end, while the offset seen from the other end
        drifts by one per nucleotide of read length. The end with the tighter
        spread is the one the protocol pins down, and it is the end to use with a
        caller that applies a single global offset.
        """
        offsets = [s.offset for s in self.scores if s.usable and s.offset is not None]
        if len(offsets) < 2:
            return None
        return max(offsets) - min(offsets)

    @property
    def quality(self) -> float:
        """A single number summarising this end, for reporting only.

        The choice between ends is made by `prefer_read_end`, not by comparing
        this value: a weighted sum saturates in high-signal data and then lets an
        irrelevant term decide.
        """
        if not self.recommendation.has_recommendation:
            return 0.0
        sharpness = min(np.log10(max(self.recommendation.sharpness, 1.0)) / 2.0, 1.0)
        frame_zero = _supported_frame_zero(
            self.recommendation.frame_fractions, self.recommendation.frame_reads
        )
        frame = max(0.0, (frame_zero - 1 / 3) / (1 - 1 / 3))
        covered = min(self.recommendation.covered_fraction / 0.5, 1.0)
        return float(0.6 * sharpness + 0.25 * frame + 0.15 * covered)


def prefer_read_end(a: "EndComparison", b: "EndComparison") -> int:
    """Order two read ends best first, as a cmp function.

    Deliberately lexicographic rather than a weighted sum, and led by offset
    spread rather than by sharpness. See the comment in the body for why
    sharpness cannot carry this decision.
    """
    if a.recommendation.has_recommendation != b.recommendation.has_recommendation:
        return -1 if a.recommendation.has_recommendation else 1
    if not a.recommendation.has_recommendation:
        return 0

    # Offset spread leads, because it is the only criterion that reflects a real
    # difference between the ends. Per read length the 5' and 3' profiles are
    # shifted copies of each other and therefore equally sharp, so a difference
    # in pooled sharpness mostly records how many read lengths each end's greedy
    # search happened to pool, not which end is better.
    spread_a, spread_b = a.offset_spread, b.offset_spread
    if spread_a is not None and spread_b is not None and spread_a != spread_b:
        return -1 if spread_a < spread_b else 1

    frame_a = _supported_frame_zero(a.recommendation.frame_fractions, a.recommendation.frame_reads)
    frame_b = _supported_frame_zero(b.recommendation.frame_fractions, b.recommendation.frame_reads)
    if abs(frame_a - frame_b) > FRAME_TIE_MARGIN:
        return -1 if frame_a > frame_b else 1

    covered_a = a.recommendation.covered_fraction
    covered_b = b.recommendation.covered_fraction
    if covered_a != covered_b:
        return -1 if covered_a > covered_b else 1

    # Last resort only, for the reason given above.
    sharp_a, sharp_b = a.recommendation.sharpness, b.recommendation.sharpness
    if max(sharp_a, sharp_b) > min(sharp_a, sharp_b) * SHARPNESS_TIE_FACTOR:
        return -1 if sharp_a > sharp_b else 1
    return 0


def compare_read_ends(
    profiles_by_end: dict[str, dict[int, np.ndarray]],
    coordinates: np.ndarray,
    read_totals: dict[int, int] | None = None,
    anchor: str = "start",
) -> tuple[EndComparison | None, list[EndComparison]]:
    """Score every read end and choose one for the selected boundary signal.

    Which end is sharper is organism and protocol dependent, so both are analysed
    and the loser is reported alongside the winner rather than discarded: seeing
    that one end is dramatically better is itself informative about the library.

    Returns (best, all_comparisons). `best` is None when no end yields a usable
    recommendation.
    """
    _validate_anchor(anchor)
    comparisons = []
    for name, profiles in profiles_by_end.items():
        read_end = read_end_for_anchor(name, anchor)
        if not profiles:
            comparisons.append(
                EndComparison(name, [], _no_recommendation([], read_end, anchor))
            )
            continue
        scores = score_read_lengths(profiles, coordinates, read_totals, read_end, anchor)
        recommendation = recommend_read_lengths(scores, profiles, coordinates, read_end, anchor)
        comparisons.append(EndComparison(name, scores, recommendation))

    comparisons.sort(key=functools.cmp_to_key(prefer_read_end))
    best = comparisons[0] if comparisons and comparisons[0].recommendation.has_recommendation else None
    return best, comparisons


def describe_end_choice(
    best: EndComparison | None, comparisons: list[EndComparison], anchor: str = "start"
) -> list[str]:
    """Explain, in words, why one read end was preferred over the other."""
    signal = _signal_for_anchor(anchor)
    if best is None:
        return [f"Neither read end produced a usable {signal} signal."]

    spread = best.offset_spread
    if spread is not None:
        lines = [
            f"{best.read_end} mapping was chosen: its estimated offset is consistent "
            f"across read lengths (varying by {spread} nt), which is the signature of "
            "the end this protocol defines precisely."
        ]
    else:
        lines = [
            f"{best.read_end} mapping was chosen; too few read lengths are usable to "
            "compare the two ends on offset consistency."
        ]

    for other in comparisons:
        if other.read_end == best.read_end:
            continue
        if not other.recommendation.has_recommendation:
            lines.append(
                f"{other.read_end} mapping produced no usable read length, so it was not used."
            )
            continue

        other_spread = other.offset_spread
        if spread is not None and other_spread is not None and other_spread > spread:
            lines.append(
                f"Seen from the {other.read_end} end the offset drifts by {other_spread} nt "
                "across read lengths, as expected when that end is the ragged one. Both "
                "ends give equally sharp profiles per read length, since for a read of "
                "fixed length they are the same profile shifted, so offset consistency "
                "rather than peak height distinguishes them."
            )
        else:
            lines.append(
                f"{other.read_end} mapping is an equally defensible choice here "
                f"({other.recommendation.sharpness:.1f}x background, "
                f"{other.recommendation.frame_bias:.0%} in the dominant frame)."
            )
    return lines
