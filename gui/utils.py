import os
import sys
import platform
from typing import Tuple


def resource_path(relative_path: str) -> str:
    """
    Obtains the absolute path of a file or resource relative to the main executable's path.

    This is particularly useful when packaging applications with PyInstaller, where resources
    are extracted to a temporary directory.

    Parameters:
        relative_path (str): The relative path to the resource.

    Returns:
        str: The absolute path to the resource.
    """
    try:
        # sys._MEIPASS is set by PyInstaller to the path of the temporary directory
        # where the application's assets are extracted.
        base_path = sys._MEIPASS  # type: ignore[attr-defined]
    except AttributeError:
        base_path = os.path.abspath(".")

    return os.path.join(base_path, relative_path)


def get_font() -> Tuple[str, int]:
    """
    Determines the default font and size based on the operating system.

    The function selects a monospace font suitable for the current operating system:
    - Windows: Consolas
    - macOS (Darwin): Menlo
    - Linux: Mono
    - Other systems: Courier (default fallback)

    Returns:
        Tuple[str, int]: A tuple containing the font name and font size.
    """
    system = platform.system()

    if system == "Windows":
        return "Consolas", 12

    elif system == "Darwin":
        return "Menlo", 12

    elif system == "Linux":
        return "Mono", 12

    else:
        return "Courier", 12