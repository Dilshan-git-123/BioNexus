"""
BioNexus FASTA Parser Module
============================
Provides robust, high-performance parsing of FASTA files formatted in
standard biological representations (.fa, .fasta, .fna, .txt).

Supports:
- Multi-line sequence aggregation
- Multiple sequence records per file
- Header decomposition (ID vs description)
- Tolerant whitespace and newline handling (CRLF, LF, CR)
- Empty record detection
- Error reporting for non-FASTA formats
"""

from dataclasses import dataclass
from io import BytesIO, StringIO, TextIOBase
from pathlib import Path
from typing import List, Union


@dataclass
class FastaRecord:
    """Represents a single FASTA sequence record."""
    id: str
    description: str
    sequence: str
    raw_header: str

    @property
    def length(self) -> int:
        """Returns the length of the nucleotide sequence."""
        return len(self.sequence)

    def __repr__(self) -> str:
        return f"FastaRecord(id='{self.id}', length={self.length})"


class FastaParseError(ValueError):
    """Exception raised when a FASTA file is malformed or invalid."""
    pass


def parse_fasta(
    source: Union[str, bytes, StringIO, BytesIO, TextIOBase, Path]
) -> List[FastaRecord]:
    """
    Parses FASTA data from a variety of sources into a list of FastaRecord objects.

    Parameters
    ----------
    source : Union[str, bytes, StringIO, BytesIO, TextIOBase, Path]
        Can be:
        - A string representing the file path or raw FASTA content
        - A bytes object (e.g. from file upload)
        - A file-like object (StringIO, BytesIO, Streamlit UploadedFile)
        - A pathlib.Path object

    Returns
    -------
    List[FastaRecord]
        List of parsed sequence records.

    Raises
    ------
    FastaParseError
        If the file content does not follow FASTA conventions or is empty.
    """
    lines: List[str] = []

    # 1. Resolve source to a list of text lines
    if isinstance(source, Path):
        with open(source, "r", encoding="utf-8", errors="replace") as f:
            lines = f.read().splitlines()
    elif isinstance(source, bytes):
        text = source.decode("utf-8", errors="replace")
        lines = text.splitlines()
    elif isinstance(source, (StringIO, BytesIO, TextIOBase)):
        content = source.read()
        if isinstance(content, bytes):
            content = content.decode("utf-8", errors="replace")
        lines = content.splitlines()
    elif isinstance(source, str):
        # Determine if string is an existing file path or raw FASTA content
        p = Path(source)
        try:
            if p.is_file():
                with open(p, "r", encoding="utf-8", errors="replace") as f:
                    lines = f.read().splitlines()
            else:
                lines = source.splitlines()
        except OSError:
            lines = source.splitlines()
    elif hasattr(source, "getvalue"):  # e.g., Streamlit UploadedFile
        raw_val = source.getvalue()
        if isinstance(raw_val, bytes):
            text = raw_val.decode("utf-8", errors="replace")
        else:
            text = str(raw_val)
        lines = text.splitlines()
    elif hasattr(source, "read"):
        content = source.read()
        if isinstance(content, bytes):
            content = content.decode("utf-8", errors="replace")
        lines = content.splitlines()
    else:
        raise FastaParseError(f"Unsupported source type: {type(source).__name__}")

    if not lines or not any(line.strip() for line in lines):
        raise FastaParseError("FASTA input is empty.")

    records: List[FastaRecord] = []
    current_header: str = ""
    current_seq_parts: List[str] = []
    has_seen_header: bool = False

    for line_idx, raw_line in enumerate(lines):
        line = raw_line.strip()

        # Skip empty lines
        if not line:
            continue

        # Skip legacy comment lines beginning with ';'
        if line.startswith(";"):
            continue

        if line.startswith(">"):
            # If we were previously buffering a record, save it
            if has_seen_header:
                seq = "".join(current_seq_parts).upper()
                rec = _build_record(current_header, seq)
                records.append(rec)
                current_seq_parts = []

            has_seen_header = True
            current_header = line[1:].strip()
        else:
            if not has_seen_header:
                # First non-empty line was not a header
                raise FastaParseError(
                    f"Invalid FASTA format: First non-blank line must begin with '>', "
                    f"found line {line_idx + 1}: '{line[:30]}...'"
                )
            # Remove all whitespace inside sequence lines
            cleaned_seq_chunk = "".join(line.split())
            current_seq_parts.append(cleaned_seq_chunk)

    # Save the final record
    if has_seen_header:
        seq = "".join(current_seq_parts).upper()
        rec = _build_record(current_header, seq)
        records.append(rec)
    else:
        raise FastaParseError("No valid sequence headers ('>') found in input.")

    return records


def _build_record(raw_header: str, sequence: str) -> FastaRecord:
    """Splits raw header into ID and description, and creates a FastaRecord."""
    if not raw_header:
        rec_id = "unnamed_sequence"
        description = ""
    else:
        tokens = raw_header.split(maxsplit=1)
        rec_id = tokens[0]
        description = tokens[1] if len(tokens) > 1 else ""

    return FastaRecord(
        id=rec_id,
        description=description,
        sequence=sequence,
        raw_header=raw_header
    )

