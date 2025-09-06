from __future__ import annotations
from pathlib import Path
import re
from typing import Optional, Tuple
from google.oauth2 import service_account
from googleapiclient.discovery import build
from googleapiclient.http import MediaFileUpload
import streamlit as st
_SCOPES = ["https://www.googleapis.com/auth/drive.file"]
def folder_id_from_link(link: str) -> Optional[str]:
    if not link: return None
    m = re.search(r"/folders/([A-Za-z0-9_-]{10,})", link)
    if m: return m.group(1)
    m = re.search(r"[?&]id=([A-Za-z0-9_-]{10,})", link)
    if m: return m.group(1)
    return None
def _creds_from_secrets():
    data = st.secrets.get("gdrive_service_account", None)
    if data is None: raise RuntimeError("Missing secret 'gdrive_service_account'.")
    if isinstance(data, str):
        import json; data = json.loads(data)
    return service_account.Credentials.from_service_account_info(data, scopes=_SCOPES)
def get_drive_service():
    return build("drive", "v3", credentials=_creds_from_secrets(), cache_discovery=False)
def upload_file(path: Path, folder_id: str, mimetype: str|None=None, name: str|None=None) -> tuple[str,str]:
    path = Path(path)
    if not path.exists(): raise FileNotFoundError(str(path))
    svc = get_drive_service()
    metadata = {"name": name or path.name, "parents": [folder_id]}
    media = MediaFileUpload(str(path), mimetype=mimetype, resumable=True)
    f = svc.files().create(body=metadata, media_body=media, fields="id, webViewLink, webContentLink").execute()
    return f["id"], f.get("webViewLink") or f.get("webContentLink", "")
