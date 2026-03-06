"""Generate a bash/zsh completion script for dip-c."""

from dip_c.data import list_bundled_data_files

SUBCOMMANDS = [
    "seg", "con", "dedup", "reg", "bincon", "clean", "impute", "impute3",
    "clean3", "reg3", "mkcon", "color", "color2", "mgcolor", "vis", "info",
    "ard", "cv", "align", "dist", "pd", "rg", "con3", "pos", "tad", "exp",
    "force", "data-path",
]

# Flags that accept a bundled data file as their argument
DATA_FILE_FLAGS = {"-c", "-n", "-l", "-L", "-m"}


def generate_completion_script():
    """Return a bash completion script string for dip-c."""
    subcommands = " ".join(SUBCOMMANDS)
    data_files = " ".join(list_bundled_data_files())
    flags = "|".join(sorted(DATA_FILE_FLAGS))

    return f'''\
_dip_c_completions() {{
    local cur prev subcmd
    cur="${{COMP_WORDS[COMP_CWORD]}}"
    prev="${{COMP_WORDS[COMP_CWORD-1]}}"

    # Find the subcommand (first non-flag argument after dip-c)
    subcmd=""
    for ((i=1; i<COMP_CWORD; i++)); do
        if [[ "${{COMP_WORDS[i]}}" != -* ]]; then
            subcmd="${{COMP_WORDS[i]}}"
            break
        fi
    done

    # Complete subcommands when none is given yet
    if [[ -z "$subcmd" ]]; then
        COMPREPLY=($(compgen -W "{subcommands}" -- "$cur"))
        return
    fi

    # After a data-file flag, complete with bundled filenames
    case "$prev" in
        {flags})
            COMPREPLY=($(compgen -W "{data_files}" -- "$cur"))
            return
            ;;
    esac

    # Default: fall back to file completion
    COMPREPLY=($(compgen -f -- "$cur"))
}}

complete -o default -F _dip_c_completions dip-c
'''
