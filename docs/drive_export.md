# Google Drive Export (Service Account)

1. In Google Cloud Console, create a **Service Account**, and enable the **Google Drive API** on the project.
2. Create a key (JSON). In Streamlit Cloud, go to **App -> Settings -> Secrets** and add:
```
[gdrive_service_account]
type = "service_account"
project_id = "YOUR_PROJECT"
private_key_id = "…"
private_key = "-----BEGIN PRIVATE KEY-----\n...\n-----END PRIVATE KEY-----\n"
client_email = "SERVICE-ACCOUNT@YOUR_PROJECT.iam.gserviceaccount.com"
client_id = "…"
token_uri = "https://oauth2.googleapis.com/token"
# (include all other standard fields from your JSON)
```
3. In Google Drive, create a folder and **share** it with the `client_email` above (Editor).
4. In the app (Informe page), paste the folder link and click **Export to Google Drive**.
