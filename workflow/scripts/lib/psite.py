"""
P-site offset estimation and read length scoring for TIS calling.

A translation initiation site caller such as ORFBounder needs to know which read
lengths carry a usable initiation signal, and by how much to shift each of them
so that a read is attributed to the codon its ribosome actually occupied. This
module derives both from the metagene profiles the workflow already computes.

The central quantity is the offset: the distance from the mapped end of a read
to the P-site of the ribosome that produced it. With 5' mapping, a ribosome
whose P-site sits on the start codon puts the read's 5' end `offset` nucleotides
upstream of it, so the offset is read off the position of the start-codon peak.

A deliberate caveat runs through the scoring: three-nucleotide periodicity is
frequently weak or absent in bacterial Ribo-seq, far more so than in eukaryotes.
Periodicity is therefore reported and scored, but a read length is not rejected
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

    sign             +1 when the mapped end lies upstream of the P-site (5'),
                     -1 when it lies downstream (3')
    search_window    where the initiation peak is looked for, relative to the
                     start codon. A 5' end peaks upstream of it, a 3' end
                     downstream, roughly a read length away
    plausible        offsets outside this range indicate noise rather than a
                     genuine initiation peak
    """

    name: str
    sign: int
    search_window: tuple[int, int]
    plausible: tuple[int, int]

    def offset_from_peak(self, peak_position: int) -> int:
        """Convert the peak's position into a positive distance to the P-site."""
        return -self.sign * peak_position

    def psite_shift(self, offset: int) -> int:
        """How far to move a profile so that it is indexed by P-site position."""
        return self.sign * offset

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

# A read length carrying less than this share of the library is too sparse to
# base an offset on, however clean its profile looks.
MIN_ABUNDANCE_FRACTION = 0.005

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

    @property
    def score(self) -> float:
        """Combined desirability of this read length for TIS calling.

        Sharpness of the initiation peak dominates, because that is the signal a
        TIS caller acts on and the one that survives in bacterial data. Frame
        bias contributes when present. Abundance breaks ties, on a log scale so
        that a very deep read length cannot outweigh a clean one.
        """
        if not self.usable:
            return 0.0
        sharpness_term = min(self.sharpness / 10.0, 1.0)
        frame_term = max(0.0, (self.frame_bias - 1 / 3) / (1 - 1 / 3))
        abundance_term = min(1.0, np.log10(1 + self.abundance * 100) / 2.0)
        return float(0.55 * sharpness_term + 0.30 * frame_term + 0.15 * abundance_term)


def _window_slice(coordinates: np.ndarray, low: int, high: int) -> np.ndarray:
    return (coordinates >= low) & (coordinates <= high)


@dataclass
class PeakEstimate:
    """The initiation peak of one profile, and how convincing it is."""

    offset: int | None
    peak_height: float
    background: float
    background_sd: float
    sharpness: float
    z_score: float

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
) -> PeakEstimate:
    """Locate the initiation peak and describe how well it stands out.

    The offset is the positive distance from the mapped read end to the start
    codon, and is None when the tallest value in the search window does not fall
    at a plausible offset.
    """
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

    background_mask = _window_slice(coordinates, *BACKGROUND_WINDOW)
    if background_mask.any():
        background_values = profile[background_mask]
        # Median for the ratio, since it is robust to a neighbouring gene
        # leaking into the upstream window; mean and sd for the z-score.
        background = float(np.median(background_values))
        background_mean = float(background_values.mean())
        background_sd = float(background_values.std())
    else:
        background = background_mean = background_sd = 0.0

    # A flat-zero background would make every peak infinitely sharp; fall back to
    # the profile's own mean so that sharpness stays comparable across lengths.
    reference = background if background > 0 else float(profile.mean())
    sharpness = float(peak_height / reference) if reference > 0 else 0.0

    # Counts are Poisson-like, so the variance is at least the mean; use that as
    # a floor to stop an unusually smooth background inflating the z-score.
    spread = max(background_sd, np.sqrt(max(background_mean, 0.0)), 1.0)
    z_score = float((peak_height - background_mean) / spread)

    offset = read_end.offset_from_peak(peak_position)
    if not read_end.is_plausible(offset):
        offset = None

    return PeakEstimate(offset, peak_height, background, background_sd, sharpness, z_score)


def frame_fractions(
    profile: np.ndarray,
    coordinates: np.ndarray,
    offset: int,
    read_end: ReadEnd = DEFAULT_READ_END,
    body_length: int = 90,
) -> tuple[float, float, float]:
    """Share of P-sites falling in each reading frame inside the ORF body.

    Frame 0 is the frame of the start codon. The first codons are skipped: the
    initiation peak itself is far larger than the elongation signal and would
    otherwise decide the answer on its own.
    """
    profile = np.asarray(profile, dtype=float)
    coordinates = np.asarray(coordinates)

    psite_positions = coordinates + read_end.psite_shift(offset)
    body = (psite_positions >= 15) & (psite_positions < 15 + body_length)
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


def periodicity_score(
    profile: np.ndarray,
    coordinates: np.ndarray,
    offset: int,
    read_end: ReadEnd = DEFAULT_READ_END,
    body_length: int = 90,
) -> float:
    """Strength of the 3-nt component of the elongation signal, in [0, 1].

    Computed as the share of the spectrum's power that sits at a period of three
    nucleotides, which is less sensitive to a single dominant codon than the
    frame fractions are.
    """
    profile = np.asarray(profile, dtype=float)
    coordinates = np.asarray(coordinates)

    psite_positions = coordinates + read_end.psite_shift(offset)
    body = (psite_positions >= 15) & (psite_positions < 15 + body_length)
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
) -> list[ReadLengthScore]:
    """Score every read length in a start-codon metagene profile."""
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

        estimate = estimate_offset(profile, coordinates, read_end)

        reasons: list[str] = []
        usable = True

        if abundance < MIN_ABUNDANCE_FRACTION:
            usable = False
            reasons.append(
                f"carries only {abundance:.2%} of the library, below the {MIN_ABUNDANCE_FRACTION:.1%} minimum"
            )

        if estimate.offset is None:
            usable = False
            reasons.append(
                "no initiation peak within "
                f"{read_end.plausible[0]}-{read_end.plausible[1]} nt of the start codon"
            )
            fractions = (0.0, 0.0, 0.0)
            periodicity = 0.0
        else:
            fractions = frame_fractions(profile, coordinates, estimate.offset, read_end)
            periodicity = periodicity_score(profile, coordinates, estimate.offset, read_end)
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
                offset=estimate.offset,
                peak_height=estimate.peak_height,
                background=estimate.background,
                sharpness=estimate.sharpness,
                z_score=estimate.z_score,
                frame_fractions=fractions,
                frame_bias=max(fractions) if any(fractions) else 0.0,
                periodicity=periodicity,
                usable=usable,
                reasons=reasons,
            )
        )

    return scores


# --------------------------------------------------------------------------
# Combining read lengths
# --------------------------------------------------------------------------


@dataclass
class Recommendation:
    """A suggested TIS caller setup, and how much to trust it."""

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

    @property
    def has_recommendation(self) -> bool:
        return bool(self.read_lengths)

    def offset_table(self) -> list[tuple[int, int]]:
        return [(length, self.offsets[length]) for length in sorted(self.read_lengths)]


def shift_profile(profile: np.ndarray, offset: int) -> np.ndarray:
    """Shift a profile so that it is indexed by P-site rather than read end.

    A read end at coordinate c belongs to a ribosome whose P-site is at c +
    offset, so the profile moves `offset` positions to the right, zero filled.
    """
    profile = np.asarray(profile, dtype=float)
    if offset == 0:
        return profile.copy()
    shifted = np.zeros_like(profile)
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
            profiles[read_length], read_end.psite_shift(offsets[read_length])
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
) -> tuple[float, tuple[float, float, float], float]:
    """Sharpness, frame fractions and periodicity of a pooled set of lengths."""
    pooled = pool_profiles(profiles, offsets, read_lengths, read_end)

    # After offset correction the initiation peak sits at coordinate 0.
    at_start = pooled[coordinates == 0]
    peak = float(at_start[0]) if at_start.size else 0.0

    background_mask = _window_slice(coordinates, *BACKGROUND_WINDOW)
    background = float(np.median(pooled[background_mask])) if background_mask.any() else 0.0
    reference = background if background > 0 else float(pooled.mean())
    sharpness = float(peak / reference) if reference > 0 else 0.0

    # The pooled profile is already indexed by P-site, so no further shift.
    fractions = frame_fractions(pooled, coordinates, 0, read_end)
    periodicity = periodicity_score(pooled, coordinates, 0, read_end)
    return sharpness, fractions, periodicity


def recommend_read_lengths(
    scores: list[ReadLengthScore],
    start_profiles: dict[int, np.ndarray],
    coordinates: np.ndarray,
    read_end: ReadEnd = DEFAULT_READ_END,
) -> Recommendation:
    """Choose the read lengths and offsets a TIS caller should be given.

    Read lengths are added greedily, best first, and a length is kept only if it
    improves the pooled initiation peak. That is preferable to taking every
    usable length: a length with a marginal peak dilutes the signal even though
    it passes on its own.
    """
    coordinates = np.asarray(coordinates)
    usable = sorted([s for s in scores if s.usable], key=lambda s: s.score, reverse=True)

    if not usable:
        return _no_recommendation(scores, read_end)

    offsets = {s.read_length: s.offset for s in usable}
    total_reads = sum(s.total_reads for s in scores) or 1

    selected = [usable[0].read_length]
    best_sharpness, best_fractions, best_periodicity = _evaluate_combination(
        start_profiles, coordinates, offsets, selected, read_end
    )

    rationale = [
        f"read length {usable[0].read_length} has the strongest initiation peak "
        f"({usable[0].sharpness:.1f}x background at offset {usable[0].offset})"
    ]

    for candidate in usable[1:]:
        trial = selected + [candidate.read_length]
        sharpness, fractions, periodicity = _evaluate_combination(
            start_profiles, coordinates, offsets, trial, read_end
        )
        if sharpness > best_sharpness * (1 + MIN_RELATIVE_GAIN):
            selected = trial
            best_sharpness, best_fractions, best_periodicity = sharpness, fractions, periodicity
            rationale.append(
                f"adding read length {candidate.read_length} (offset {candidate.offset}) "
                f"raises the pooled peak to {sharpness:.1f}x"
            )
        else:
            rationale.append(
                f"read length {candidate.read_length} was left out: it moves the pooled "
                f"peak from {best_sharpness:.1f}x to {sharpness:.1f}x, short of the "
                f"{MIN_RELATIVE_GAIN:.0%} gain required to include it"
            )

    covered = sum(s.total_reads for s in scores if s.read_length in selected) / total_reads
    confidence, warnings = _assess_confidence(best_sharpness, best_fractions, best_periodicity, covered)

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
    )


def _no_recommendation(
    scores: list[ReadLengthScore], read_end: ReadEnd = DEFAULT_READ_END
) -> Recommendation:
    """Explain why no read length is usable instead of inventing a setup."""
    warnings = [
        "No read length shows a usable translation initiation signal, so no TIS caller "
        "setup is suggested."
    ]
    if not scores:
        warnings.append("No reads were found in the metagene window at all.")
    else:
        best = max(scores, key=lambda s: s.sharpness, default=None)
        if best is not None and best.sharpness > 0:
            warnings.append(
                f"The closest candidate was read length {best.read_length} at "
                f"{best.sharpness:.1f}x background ({best.z_score:.1f} standard deviations), "
                f"below the required {MIN_PEAK_RATIO:.0f}x and {MIN_PEAK_Z:.0f}."
            )
        warnings.append(
            "Common causes: the library is RNA-seq rather than Ribo-seq, the reads are "
            "dominated by rRNA or tRNA, or the annotation start codons are inaccurate."
        )
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
    )


def _assess_confidence(
    sharpness: float,
    fractions: tuple[float, float, float],
    periodicity: float,
    covered: float,
) -> tuple[str, list[str]]:
    """Grade a recommendation, and say which evidence it rests on."""
    frame_bias = max(fractions) if any(fractions) else 0.0
    in_frame = fractions[0] if any(fractions) else 0.0
    warnings: list[str] = []

    if sharpness >= 5.0 and frame_bias >= 0.5:
        confidence = "high"
    elif sharpness >= 5.0:
        confidence = "medium"
        warnings.append(
            f"The initiation peak is clear ({sharpness:.1f}x) but the reading frame bias is "
            f"weak ({frame_bias:.0%} in the dominant frame). This is common in bacterial "
            "Ribo-seq; the offsets rest on the start-codon peak alone."
        )
    elif sharpness >= 3.0:
        confidence = "medium"
    else:
        confidence = "low"
        warnings.append(
            f"The pooled initiation peak is only {sharpness:.1f}x background. Treat the "
            "offsets as provisional and inspect the metagene plots before relying on them."
        )

    if frame_bias >= 0.5 and in_frame < frame_bias:
        dominant = int(np.argmax(fractions))
        warnings.append(
            f"The dominant reading frame is {dominant}, not 0. The estimated offsets may be "
            f"off by {dominant} nt, or the annotated start codons may be shifted."
        )

    if covered < 0.3:
        warnings.append(
            f"The selected read lengths cover only {covered:.0%} of the library. Check the "
            "read length distribution for a broader usable range."
        )

    if periodicity < 0.15 and confidence != "low":
        warnings.append(
            f"Three-nucleotide periodicity is weak ({periodicity:.0%} of spectral power). "
            "Sub-codon assignment will be unreliable even though the initiation signal is good."
        )

    return confidence, warnings


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
        frame = max(0.0, (self.recommendation.frame_bias - 1 / 3) / (1 - 1 / 3))
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

    frame_a, frame_b = a.recommendation.frame_bias, b.recommendation.frame_bias
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
) -> tuple[EndComparison | None, list[EndComparison]]:
    """Score every read end and pick the one with the better initiation signal.

    Which end is sharper is organism and protocol dependent, so both are analysed
    and the loser is reported alongside the winner rather than discarded: seeing
    that one end is dramatically better is itself informative about the library.

    Returns (best, all_comparisons). `best` is None when no end yields a usable
    recommendation.
    """
    comparisons = []
    for name, profiles in profiles_by_end.items():
        read_end = READ_ENDS[name]
        if not profiles:
            comparisons.append(
                EndComparison(name, [], _no_recommendation([], read_end))
            )
            continue
        scores = score_read_lengths(profiles, coordinates, read_totals, read_end)
        recommendation = recommend_read_lengths(scores, profiles, coordinates, read_end)
        comparisons.append(EndComparison(name, scores, recommendation))

    comparisons.sort(key=functools.cmp_to_key(prefer_read_end))
    best = comparisons[0] if comparisons and comparisons[0].recommendation.has_recommendation else None
    return best, comparisons


def describe_end_choice(best: EndComparison | None, comparisons: list[EndComparison]) -> list[str]:
    """Explain, in words, why one read end was preferred over the other."""
    if best is None:
        return ["Neither read end produced a usable initiation signal."]

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
