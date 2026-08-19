from __future__ import annotations

from pathlib import Path

from txnova.staging import STAGING_PREFIX, StagingDir, remove_orphan_staging


def test_orphan_dead_pid_removed(tmp_path: Path) -> None:
    orphan = tmp_path / f"{STAGING_PREFIX}999999999"
    orphan.mkdir()
    (orphan / "partial.txt").write_text("x")
    removed = remove_orphan_staging(tmp_path)
    assert "999999999" in removed[0]
    assert not orphan.exists()


def test_live_pid_kept(tmp_path: Path) -> None:
    import os

    live = tmp_path / f"{STAGING_PREFIX}{os.getpid()}"
    live.mkdir()
    (live / "partial.txt").write_text("x")
    remove_orphan_staging(tmp_path)
    assert live.exists()
    assert (live / "partial.txt").exists()


def test_publish_then_cleanup(tmp_path: Path) -> None:
    with StagingDir(tmp_path) as st:
        st.write_text("preflight.json", '{"ok": true}\n')
        st.publish(["preflight.json"])
    assert (tmp_path / "preflight.json").read_text() == '{"ok": true}\n'
    leftover = list(tmp_path.glob(f"{STAGING_PREFIX}*"))
    assert leftover == []


def test_failure_does_not_publish(tmp_path: Path) -> None:
    try:
        with StagingDir(tmp_path) as st:
            st.write_text("preflight.json", "partial")
            raise RuntimeError("boom")
    except RuntimeError:
        pass
    assert not (tmp_path / "preflight.json").exists()
    leftover = list(tmp_path.glob(f"{STAGING_PREFIX}*"))
    assert leftover == []
