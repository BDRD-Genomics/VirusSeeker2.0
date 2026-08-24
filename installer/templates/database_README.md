# VirusSeeker database mount

Database profile is recorded in `../config/run.env`.

Docker mode mounts this entire directory read-only:

```bash
-v <VS_ROOT>/databases:/database:ro
```

Stable paths used by `config/VS.config`:

| Host path | Container path |
|---|---|
| `ref/GRCh38.p14` | `/database/ref/GRCh38.p14` |
| `nr/nr` | `/database/nr/nr` |
| `nt/nt` | `/database/nt/nt` |
| `core_nt/core_nt` | `/database/core_nt/core_nt` |
| `nr_dmnd/nr.dmnd` | `/database/nr_dmnd/nr.dmnd` |
| `core_nt_mmseqs/core_nt` | `/database/core_nt_mmseqs/core_nt` |
| `VirusDBNR/virus_nr.dmnd` | `/database/VirusDBNR/virus_nr.dmnd` |
| `VirusDBNT/virus_nt` | `/database/VirusDBNT/virus_nt` |
| `taxdump/` | `/database/taxdump/` |
| `vhunter_acc.db` | `/database/vhunter_acc.db` |
| `taxa.sqlite` | `/database/taxa.sqlite` |
| `famdb/*.h5` | `/database/famdb/*.h5` |

Dated database releases are stored under `releases/`. Stable symlinks point to
the active release, allowing updates without changing Docker mounts or
`VS.config`.
