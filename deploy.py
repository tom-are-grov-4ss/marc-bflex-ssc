"""Deploy Python Environment"""

import venv
from pathlib import Path
from subprocess import run

CURRENT_DIR = Path(__file__).parent.resolve()
BASENAME = CURRENT_DIR.name

venv.create(CURRENT_DIR / ".venv", prompt=BASENAME, with_pip=True)

run(
    [
        str(CURRENT_DIR / ".venv" / "Scripts" / "python"),
        "-m",
        "pip",
        "install",
        "--upgrade",
        "pip",
    ]
)

run(
    [
        str(CURRENT_DIR / ".venv" / "Scripts" / "python"),
        "-m",
        "pip",
        "install",
        "-r",
        str(CURRENT_DIR / "requirements.txt"),
    ]
)
