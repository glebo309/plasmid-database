from flask import Flask, render_template_string, request, jsonify, send_from_directory, session, redirect, url_for, Response
import sqlite3
import os
import shutil
import json
import math
import re
import io
import zipfile
import tempfile
from datetime import datetime
from threading import Thread
import time
import schedule
from werkzeug.utils import secure_filename
from werkzeug.security import generate_password_hash, check_password_hash
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

app = Flask(__name__)
app.config['MAX_CONTENT_LENGTH'] = 50 * 1024 * 1024  # 50MB max file size
app.secret_key = os.urandom(24)  # Generate random secret key for sessions

# Allow iframe embedding
@app.after_request
def after_request(response):
    response.headers.pop('X-Frame-Options', None)
    response.headers['Content-Security-Policy'] = "frame-ancestors *"
    response.headers['Cache-Control'] = 'no-store'
    return response

# Admin password (change this!)
ADMIN_PASSWORD_HASH = generate_password_hash('admin123')  # Default password: admin123

# Database file location
DB_FILE = 'plasmid_library.db'
BACKUP_DIR = 'backups'
UPLOAD_DIR = 'uploads'

# Ensure directories exist
os.makedirs(BACKUP_DIR, exist_ok=True)
os.makedirs(UPLOAD_DIR, exist_ok=True)

ALLOWED_EXTENSIONS = {'txt', 'pdf', 'png', 'jpg', 'jpeg', 'gif', 'dna', 'seq', 'gb', 'ape', 'fasta', 'docx', 'xlsx'}

def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

# Initialize database
def init_db():
    conn = sqlite3.connect(DB_FILE)
    c = conn.cursor()
    c.execute('''CREATE TABLE IF NOT EXISTS plasmids
                 (id INTEGER PRIMARY KEY AUTOINCREMENT,
                  name TEXT NOT NULL,
                  number INTEGER UNIQUE,
                  plasmid_backbone TEXT,
                  inductor TEXT,
                  antibiotic_resistance TEXT,
                  reference TEXT,
                  comment TEXT,
                  created_at TEXT DEFAULT CURRENT_TIMESTAMP,
                  modified_at TEXT DEFAULT CURRENT_TIMESTAMP,
                  files TEXT,
                  added_on TEXT,
                  added_by TEXT)''')
    
    # Add columns missing from older databases
    for column in ('files', 'added_on', 'added_by'):
        try:
            c.execute(f'ALTER TABLE plasmids ADD COLUMN {column} TEXT')
            conn.commit()
        except sqlite3.OperationalError:
            pass  # Column already exists
    
    # Audit trail table - tracks ALL changes
    c.execute('''CREATE TABLE IF NOT EXISTS audit_log
                 (id INTEGER PRIMARY KEY AUTOINCREMENT,
                  timestamp TEXT NOT NULL,
                  action TEXT NOT NULL,
                  record_id INTEGER,
                  old_data TEXT,
                  new_data TEXT,
                  ip_address TEXT)''')
    conn.commit()
    conn.close()

def log_audit(action, record_id=None, old_data=None, new_data=None, ip_address=None):
    """Log all database changes to audit trail"""
    conn = sqlite3.connect(DB_FILE)
    c = conn.cursor()
    c.execute('''INSERT INTO audit_log (timestamp, action, record_id, old_data, new_data, ip_address)
                 VALUES (?, ?, ?, ?, ?, ?)''',
              (datetime.now().isoformat(), action, record_id, 
               str(old_data) if old_data else None,
               str(new_data) if new_data else None,
               ip_address))
    conn.commit()
    conn.close()
    print(f"[AUDIT] {action} - Record ID: {record_id} - IP: {ip_address}")

def backup_database():
    """Create timestamped backup of database"""
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    backup_file = os.path.join(BACKUP_DIR, f'plasmid_library_backup_{timestamp}.db')
    try:
        shutil.copy2(DB_FILE, backup_file)
        print(f"[BACKUP] Database backed up to: {backup_file}")
        log_audit('BACKUP_CREATED', new_data=backup_file, ip_address='SYSTEM')
        
        # Keep only last 12 backups (3 months of Fridays)
        backups = sorted([f for f in os.listdir(BACKUP_DIR) if f.endswith('.db')])
        if len(backups) > 12:
            for old_backup in backups[:-12]:
                os.remove(os.path.join(BACKUP_DIR, old_backup))
                print(f"[BACKUP] Deleted old backup: {old_backup}")
    except Exception as e:
        print(f"[ERROR] Backup failed: {e}")

def schedule_backups():
    """Schedule automatic backups every Friday at 6 PM"""
    schedule.every().friday.at("18:00").do(backup_database)
    
    while True:
        schedule.run_pending()
        time.sleep(3600)  # Check every hour

# HTML Template for main view (NO DELETE, only EDIT)
HTML_TEMPLATE = '''
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Plasmid Library</title>
<style>
:root { --fg: #1a1a1a; --fg2: #555; --fg3: #999; --bg: #fff; --bg2: #fafafa; --border: #e0e0e0; --accent: #2e7d4f; --accent-light: #eaf6ef; --danger: #c0392b; --mono: "SF Mono", "Cascadia Code", "Fira Code", Consolas, monospace; }
*, *::before, *::after { margin: 0; padding: 0; box-sizing: border-box; }
html { font-size: 14px; }
body { font-family: "Helvetica Neue", Helvetica, Arial, sans-serif; color: var(--fg); background: var(--bg); line-height: 1.5; }

/* top bar */
.topbar { display: flex; align-items: center; justify-content: space-between; padding: 8px 20px; border-bottom: 1px solid var(--border); background: var(--bg); position: sticky; top: 0; z-index: 100; }
.topbar h1 { font-size: 16px; font-weight: 600; letter-spacing: -0.01em; }
.topbar-actions { display: flex; gap: 6px; align-items: center; }

/* buttons */
.btn { display: inline-flex; align-items: center; gap: 4px; padding: 6px 14px; font-size: 12px; font-weight: 500; border: 1px solid var(--border); border-radius: 3px; background: var(--bg); color: var(--fg); cursor: pointer; text-decoration: none; white-space: nowrap; transition: background 0.1s; }
.btn:hover { background: var(--bg2); }
.btn-primary { background: var(--fg); color: var(--bg); border-color: var(--fg); }
.btn-primary:hover { background: #333; }
.btn-ghost { border-color: transparent; }
.btn-ghost:hover { background: var(--bg2); }
.btn-sm { padding: 3px 8px; font-size: 11px; }
.btn-danger { color: var(--danger); border-color: var(--danger); }
.btn-danger:hover { background: #fdf0ef; }

/* content area */
.content { padding: 10px 20px; }

/* toolbar */
.toolbar { display: flex; gap: 6px; margin-bottom: 8px; flex-wrap: wrap; align-items: center; }
.search-input { flex: 1; min-width: 240px; padding: 6px 10px; font-size: 12px; border: 1px solid var(--border); border-radius: 3px; background: var(--bg); font-family: inherit; }
.search-input:focus { outline: none; border-color: var(--accent); }
.filter-select { padding: 6px 8px; font-size: 12px; border: 1px solid var(--border); border-radius: 3px; background: var(--bg); font-family: inherit; cursor: pointer; }
.filter-select:focus { outline: none; border-color: var(--accent); }
.stats { font-size: 10px; color: var(--fg3); margin-bottom: 6px; }

/* add form */
.form-panel { display: none; padding: 16px; margin-bottom: 16px; border: 1px solid var(--border); background: var(--bg2); }
.form-panel.active { display: block; }
.form-panel h3 { font-size: 13px; font-weight: 600; margin-bottom: 12px; }
.form-grid { display: grid; grid-template-columns: repeat(3, 1fr); gap: 12px; }
.form-group { margin-bottom: 8px; }
.form-group.full-width { grid-column: 1 / -1; }
.form-group label { display: block; font-size: 11px; font-weight: 600; color: var(--fg2); margin-bottom: 4px; text-transform: uppercase; letter-spacing: 0.03em; }
.form-group input, .form-group textarea { width: 100%; padding: 6px 8px; font-size: 12px; border: 1px solid var(--border); border-radius: 3px; font-family: inherit; background: var(--bg); }
.form-group input:focus, .form-group textarea:focus { outline: none; border-color: var(--accent); }
.form-group textarea { min-height: 60px; resize: vertical; }
.form-group input[readonly] { background: var(--bg2); color: var(--fg3); }
.form-group small { font-size: 11px; color: var(--fg3); }
.form-actions { display: flex; gap: 6px; margin-top: 8px; }

/* table */
.table-wrap { }
table { width: 100%; border-collapse: collapse; }
th { padding: 5px 8px; text-align: left; font-size: 12px; font-weight: 600; color: var(--fg2); text-transform: uppercase; letter-spacing: 0.03em; border-bottom: 2px solid var(--fg); cursor: pointer; user-select: none; white-space: nowrap; background: var(--bg); position: sticky; top: 42px; z-index: 10; }
th:hover { color: var(--fg); }
th.sort-asc::after { content: " \\2191"; }
th.sort-desc::after { content: " \\2193"; }
td { padding: 5px 8px; border-bottom: 1px solid var(--border); border-right: 1px solid #eee; font-size: 13px; vertical-align: top; }
td:last-child { border-right: none; }
tbody tr:hover { background: var(--accent-light); }
tbody tr.expanded .col-narrow, tbody tr.expanded .col-comment, tbody tr.expanded .col-minor { white-space: normal; text-overflow: unset; word-break: break-word; }
tbody tr.hidden { display: none; }
.col-id { width: 40px; font-weight: 600; font-family: var(--mono); font-size: 10px; max-width: 40px; }
.col-name { min-width: 120px; max-width: 220px; font-weight: 500; }
.col-narrow { max-width: 80px; overflow: hidden; text-overflow: ellipsis; white-space: nowrap; }
.col-comment { max-width: 160px; overflow: hidden; text-overflow: ellipsis; white-space: nowrap; }
.col-minor { max-width: 100px; overflow: hidden; text-overflow: ellipsis; white-space: nowrap; font-size: 11px; color: var(--fg2); }
.col-action { width: 50px; text-align: center; max-width: 50px; }
td.muted { color: var(--fg3); }

/* file chips */
td:has(.file-list) { white-space: normal; max-width: 200px; }
.file-list { display: flex; flex-direction: column; gap: 2px; }
.file-chip { display: inline-flex; align-items: center; gap: 4px; font-size: 10px; }
.file-chip a { color: var(--accent); text-decoration: none; }
.file-chip a:hover { text-decoration: underline; }
.file-chip .actions { display: inline-flex; gap: 2px; }
.file-chip .actions button, .file-chip .actions a { font-size: 10px; padding: 1px 6px; border: 1px solid var(--border); border-radius: 2px; background: var(--bg); color: var(--fg2); cursor: pointer; text-decoration: none; display: inline-block; }
.file-chip .actions button:hover, .file-chip .actions a:hover { background: var(--bg2); }

/* modal */
.modal { display: none; position: fixed; z-index: 1000; inset: 0; background: rgba(0,0,0,0.4); }
.modal.active { display: flex; align-items: center; justify-content: center; }
.modal-box { background: var(--bg); border: 1px solid var(--border); width: 90%; max-width: 720px; max-height: 90vh; overflow-y: auto; }
.modal-head { display: flex; justify-content: space-between; align-items: center; padding: 12px 16px; border-bottom: 1px solid var(--border); }
.modal-head h2 { font-size: 14px; font-weight: 600; }
.modal-body { padding: 16px; }
.modal-foot { display: flex; gap: 6px; padding: 12px 16px; border-top: 1px solid var(--border); }

/* file items in edit modal */
.file-tag { display: inline-flex; align-items: center; gap: 6px; padding: 4px 8px; margin: 2px; background: var(--bg2); border: 1px solid var(--border); font-size: 11px; }
.file-tag a { color: var(--accent); text-decoration: none; }
.file-tag .remove-btn { background: none; border: none; color: var(--fg3); cursor: pointer; font-size: 14px; line-height: 1; padding: 0 2px; }
.file-tag .remove-btn:hover { color: var(--danger); }

/* sequence viewer modal */
#seqModal .modal-box { max-width: 1000px; }
.seq-info { display: flex; gap: 12px; margin-bottom: 12px; flex-wrap: wrap; }
.seq-badge { padding: 4px 10px; background: var(--bg2); border: 1px solid var(--border); font-size: 11px; }
.seq-badge.highlight { background: #fdf0ef; border-color: var(--danger); }
.seq-tabs { display: flex; border-bottom: 1px solid var(--border); margin-bottom: 12px; }
.seq-tab { padding: 8px 16px; font-size: 12px; border: none; background: none; cursor: pointer; color: var(--fg2); border-bottom: 2px solid transparent; }
.seq-tab:hover { color: var(--fg); }
.seq-tab.active { color: var(--accent); border-bottom-color: var(--accent); font-weight: 500; }
.seq-map { text-align: center; min-height: 300px; }
.seq-text-box { max-height: 400px; overflow-y: auto; background: #1c1c1c; padding: 12px; }
.seq-text { font-family: var(--mono); font-size: 11px; color: #ccc; word-break: break-all; line-height: 1.6; white-space: pre-wrap; }
.seq-text .insert-highlight { background: rgba(192,57,43,0.25); color: #e74c3c; font-weight: bold; }
.seq-text .line-num { color: #666; margin-right: 8px; user-select: none; }
.feature-table { width: 100%; border-collapse: collapse; font-size: 11px; }
.feature-table th { padding: 6px 8px; text-align: left; font-size: 10px; font-weight: 600; color: var(--fg2); text-transform: uppercase; border-bottom: 2px solid var(--border); background: var(--bg); }
.feature-table td { padding: 6px 8px; border-bottom: 1px solid #f0f0f0; }
.feature-type { display: inline-block; padding: 1px 6px; font-size: 10px; font-weight: 600; color: #fff; border-radius: 2px; }
.empty-state { text-align: center; padding: 60px; color: var(--fg3); }

@media (max-width: 900px) {
    .form-grid { grid-template-columns: 1fr; }
    .toolbar { flex-direction: column; }
    .search-input { min-width: 100%; }
}
</style>
</head>
<body>
<div class="topbar">
    <h1>Plasmid Library</h1>
    <div class="topbar-actions">
        <a href="/export/html" target="_blank" class="btn">Export</a>
        <a href="/admin" class="btn">Admin</a>
        <button class="btn btn-primary" onclick="toggleAddForm()">+ Add Entry</button>
    </div>
</div>

<div class="content">
    <div class="form-panel" id="addForm">
        <h3>Add New Plasmid</h3>
        <form onsubmit="addPlasmid(event)" enctype="multipart/form-data">
            <div class="form-grid">
                <div class="form-group"><label>ID</label><input type="number" id="nextIdField" name="number" readonly required></div>
                <div class="form-group"><label>Name *</label><input type="text" name="name" required></div>
                <div class="form-group"><label>Plasmid Backbone</label><input type="text" name="plasmid_backbone"></div>
                <div class="form-group"><label>Inductor</label><input type="text" name="inductor"></div>
                <div class="form-group"><label>Antibiotic Resistance</label><input type="text" name="antibiotic_resistance"></div>
                <div class="form-group"><label>Reference</label><input type="text" name="reference"></div>
                <div class="form-group"><label>Added By *</label><input type="text" name="added_by" required placeholder="Your name"></div>
                <div class="form-group full-width"><label>Comment</label><textarea name="comment"></textarea></div>
                <div class="form-group full-width"><label>Sequence File</label><input type="file" name="files" multiple id="addFileInput"><small>Optional (.gb, .seq, .fasta, .ape, .dna)</small></div>
            </div>
            <div class="form-actions">
                <button type="submit" class="btn btn-primary">Save</button>
                <button type="button" class="btn" onclick="toggleAddForm()">Cancel</button>
            </div>
        </form>
    </div>

    <div class="toolbar">
        <input type="text" class="search-input" id="searchInput" placeholder="Search..." onkeyup="filterTable()">
        <select class="filter-select" id="backboneFilter" onchange="filterTable()"><option value="">All Backbones</option></select>
        <select class="filter-select" id="resistanceFilter" onchange="filterTable()"><option value="">All Resistances</option></select>
        <select class="filter-select" id="inductorFilter" onchange="filterTable()"><option value="">All Inductors</option></select>
    </div>

    <div class="stats"><span id="visibleCount">0</span> / <span id="totalCount">0</span> entries</div>

    {% if plasmids %}
    <div class="table-wrap">
        <table id="plasmidTable">
            <thead><tr>
                <th class="col-id" onclick="sortTable(0)">ID</th>
                <th class="col-name" onclick="sortTable(1)">Name</th>
                <th onclick="sortTable(2)">Backbone</th>
                <th onclick="sortTable(3)">Resistance</th>
                <th onclick="sortTable(4)">Inductor</th>
                <th>Files</th>
                <th class="col-comment" onclick="sortTable(6)">Comment</th>
                <th class="col-minor" onclick="sortTable(7)">Reference</th>
                <th class="col-minor" onclick="sortTable(8)">Added By</th>
                <th class="col-action"></th>
            </tr></thead>
            <tbody>
            {% for p in plasmids %}
            <tr data-id="{{ p[2] or 0 }}" data-name="{{ p[1] }}" data-backbone="{{ p[3] or '' }}" data-inductor="{{ p[4] or '' }}" data-resistance="{{ p[5] or '' }}">
                <td class="col-id">{{ p[2] or '' }}</td>
                <td class="col-name">{{ p[1] }}</td>
                <td{% if not p[3] %} class="muted"{% endif %} title="{{ p[3] or '' }}">{{ p[3] or '-' }}</td>
                <td{% if not p[5] %} class="muted"{% endif %} title="{{ p[5] or '' }}">{{ p[5] or '-' }}</td>
                <td class="col-narrow"{% if not p[4] %} style="color:var(--fg3)"{% endif %} title="{{ p[4] or '' }}">{{ p[4] or '-' }}</td>
                <td>{% if p[10] %}{% set files = p[10].split(',') %}<div class="file-list">{% for file in files %}{% set fname = file.strip() %}{% set display = fname.split('_', 3)[-1] if '_' in fname else fname %}{% set ext = fname.rsplit('.', 1)[-1].lower() if '.' in fname else '' %}<div class="file-chip"><a href="/download/{{ fname }}" download title="{{ display }}">{{ display }}</a><span class="actions">{% if ext in ['gb', 'gbk', 'ape', 'dna', 'fasta', 'seq'] %}<button data-filename="{{ fname }}" onclick="viewSequence(this.dataset.filename)">View</button>{% endif %}</span></div>{% endfor %}</div>{% endif %}</td>
                <td class="col-comment"{% if not p[7] %} style="color:var(--fg3)"{% endif %} title="{{ p[7] or '' }}">{{ p[7] or '-' }}</td>
                <td class="col-minor"{% if not p[6] %} style="color:var(--fg3)"{% endif %} title="{{ p[6] or '' }}">{{ p[6] or '-' }}</td>
                <td class="col-minor"{% if not p[12] %} style="color:var(--fg3)"{% endif %} title="{{ p[12] or '' }}">{{ p[12] or '-' }}</td>
                <td class="col-action"><button class="btn btn-sm btn-ghost" onclick='editPlasmid({{ p[0] }}, {{ p|tojson }})'>Edit</button></td>
            </tr>
            {% endfor %}
            </tbody>
        </table>
    </div>
    {% else %}
    <div class="empty-state"><p>No entries yet.</p></div>
    {% endif %}
</div>

<div id="editModal" class="modal">
    <div class="modal-box">
        <div class="modal-head"><h2>Edit Plasmid</h2><button class="btn btn-sm" onclick="closeEditModal()">Close</button></div>
        <form onsubmit="updatePlasmid(event)" enctype="multipart/form-data">
        <div class="modal-body">
            <input type="hidden" id="edit_id" name="id">
            <div class="form-grid">
                <div class="form-group"><label>ID</label><input type="number" id="edit_number" name="number" readonly required></div>
                <div class="form-group"><label>Name *</label><input type="text" id="edit_name" name="name" required></div>
                <div class="form-group"><label>Plasmid Backbone</label><input type="text" id="edit_plasmid_backbone" name="plasmid_backbone"></div>
                <div class="form-group"><label>Inductor</label><input type="text" id="edit_inductor" name="inductor"></div>
                <div class="form-group"><label>Antibiotic Resistance</label><input type="text" id="edit_antibiotic_resistance" name="antibiotic_resistance"></div>
                <div class="form-group"><label>Reference</label><input type="text" id="edit_reference" name="reference"></div>
                <div class="form-group"><label>Added By</label><input type="text" id="edit_added_by" name="added_by"></div>
                <div class="form-group full-width"><label>Comment</label><textarea id="edit_comment" name="comment"></textarea></div>
                <div class="form-group full-width"><label>Files</label><div id="current_files"></div></div>
                <div class="form-group full-width"><label>Add Files</label><input type="file" name="files" multiple></div>
            </div>
        </div>
        <div class="modal-foot">
            <button type="submit" class="btn btn-primary">Save Changes</button>
            <button type="button" class="btn" onclick="closeEditModal()">Cancel</button>
        </div>
        </form>
    </div>
</div>

<div id="seqModal" class="modal">
    <div class="modal-box">
        <div class="modal-head">
            <h2 id="seqModalTitle">Sequence Viewer</h2>
            <div style="display:flex;gap:6px"><a id="seqDownloadBtn" class="btn btn-sm" download>Download</a><button class="btn btn-sm" onclick="closeSeqModal()">Close</button></div>
        </div>
        <div class="modal-body">
            <div class="seq-info" id="seqInfo"></div>
            <div class="seq-tabs">
                <button class="seq-tab active" onclick="showSeqTab('map', event)">Map</button>
                <button class="seq-tab" onclick="showSeqTab('sequence', event)">Sequence</button>
                <button class="seq-tab" onclick="showSeqTab('features', event)">Features</button>
            </div>
            <div id="seqTabMap" class="seq-tab-content"><div id="seqMapContainer" class="seq-map"><p style="color:var(--fg3);padding:40px">Loading...</p></div></div>
            <div id="seqTabSequence" class="seq-tab-content" style="display:none"><div class="seq-text-box" id="seqTextContainer"><div class="seq-text" id="seqText">Loading...</div></div></div>
            <div id="seqTabFeatures" class="seq-tab-content" style="display:none"><div id="seqFeaturesContainer">Loading...</div></div>
        </div>
    </div>
</div>

<script>
let sortDirection = {};
window.onload = function() {
    populateFilters(); updateStats();
    document.querySelectorAll('#plasmidTable tbody tr').forEach(row => {
        row.addEventListener('click', function(e) {
            if (e.target.closest('button, a, input')) return;
            this.classList.toggle('expanded');
        });
    });
};

async function toggleAddForm() {
    const form = document.getElementById('addForm');
    const isActive = form.classList.toggle('active');
    if (isActive) {
        const r = await fetch('/api/next-id');
        const d = await r.json();
        document.getElementById('nextIdField').value = d.next_id;
    }
}

function populateFilters() {
    const rows = document.querySelectorAll('#plasmidTable tbody tr');
    const sets = { backbone: new Set(), resistance: new Set(), inductor: new Set() };
    rows.forEach(r => {
        ['backbone','resistance','inductor'].forEach(k => { const v = r.dataset[k]; if (v && v !== '-') sets[k].add(v); });
    });
    populateSelect('backboneFilter', Array.from(sets.backbone).sort());
    populateSelect('resistanceFilter', Array.from(sets.resistance).sort());
    populateSelect('inductorFilter', Array.from(sets.inductor).sort());
}
function populateSelect(id, values) {
    const sel = document.getElementById(id);
    const cur = sel.value;
    const label = id.replace('Filter','').replace(/([A-Z])/g,' $1');
    sel.innerHTML = '<option value="">All ' + label + 's</option>';
    values.forEach(v => { const o = document.createElement('option'); o.value = v; o.textContent = v; sel.appendChild(o); });
    sel.value = cur;
}

function filterTable() {
    const q = document.getElementById('searchInput').value.toLowerCase();
    const bf = document.getElementById('backboneFilter').value;
    const rf = document.getElementById('resistanceFilter').value;
    const inf = document.getElementById('inductorFilter').value;
    let vis = 0;
    document.querySelectorAll('#plasmidTable tbody tr').forEach(row => {
        const ok = row.textContent.toLowerCase().includes(q) && (!bf || row.dataset.backbone === bf) && (!rf || row.dataset.resistance === rf) && (!inf || row.dataset.inductor === inf);
        row.classList.toggle('hidden', !ok);
        if (ok) vis++;
    });
    updateStats(vis);
}
function updateStats(v) {
    const t = document.querySelectorAll('#plasmidTable tbody tr').length;
    document.getElementById('totalCount').textContent = t;
    document.getElementById('visibleCount').textContent = v !== undefined ? v : t;
}
function sortTable(col) {
    const table = document.getElementById('plasmidTable');
    const tbody = table.querySelector('tbody');
    const rows = Array.from(tbody.querySelectorAll('tr'));
    const dir = (sortDirection[col] || 'asc') === 'asc' ? 'desc' : 'asc';
    sortDirection = { [col]: dir };
    table.querySelectorAll('th').forEach((h,i) => { h.classList.remove('sort-asc','sort-desc'); if (i===col) h.classList.add('sort-'+dir); });
    rows.sort((a,b) => {
        let av = a.cells[col].textContent.trim(), bv = b.cells[col].textContent.trim();
        if (av==='-') av=''; if (bv==='-') bv='';
        const an=parseFloat(av), bn=parseFloat(bv);
        if (!isNaN(an)&&!isNaN(bn)) return dir==='asc'?an-bn:bn-an;
        return dir==='asc'?av.localeCompare(bv):bv.localeCompare(av);
    });
    rows.forEach(r => tbody.appendChild(r));
}

async function addPlasmid(e) {
    e.preventDefault();
    const r = await fetch('/api/plasmids', {method:'POST', body: new FormData(e.target)});
    if (r.ok) location.reload(); else { const d = await r.json().catch(()=>({})); alert(d.error||'Error'); }
}
function editPlasmid(id, data) {
    document.getElementById('edit_id').value = id;
    document.getElementById('edit_number').value = data[2]||'';
    document.getElementById('edit_name').value = data[1]||'';
    document.getElementById('edit_plasmid_backbone').value = data[3]||'';
    document.getElementById('edit_inductor').value = data[4]||'';
    document.getElementById('edit_antibiotic_resistance').value = data[5]||'';
    document.getElementById('edit_reference').value = data[6]||'';
    document.getElementById('edit_comment').value = data[7]||'';
    document.getElementById('edit_added_by').value = data[12]||'';
    window._editId = id;
    const fd = document.getElementById('current_files');
    if (data[10]) {
        fd.innerHTML = data[10].split(',').map(f => {
            const dn = f.includes('_')?f.split('_').slice(3).join('_'):f;
            const safe = f.replace(/[^a-zA-Z0-9]/g,'_');
            return '<span class="file-tag" id="file-'+safe+'"><a href="/download/'+f+'" download>'+dn+'</a><button type="button" class="remove-btn" data-fn="'+f+'">&times;</button></span>';
        }).join('');
        fd.querySelectorAll('.remove-btn').forEach(b => b.onclick = function(){ removeFile(this.dataset.fn); });
    } else { fd.innerHTML = '<span style="color:var(--fg3);font-size:11px">No files attached</span>'; }
    document.getElementById('editModal').classList.add('active');
}
function closeEditModal() { document.getElementById('editModal').classList.remove('active'); }
async function removeFile(fn) {
    if (!confirm('Remove this file?')) return;
    const r = await fetch('/api/plasmids/'+window._editId+'/remove-file', {method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({filename:fn})});
    if (r.ok) {
        const el = document.getElementById('file-'+fn.replace(/[^a-zA-Z0-9]/g,'_'));
        if (el) el.remove();
        const fd = document.getElementById('current_files');
        if (!fd.querySelector('.file-tag')) fd.innerHTML = '<span style="color:var(--fg3);font-size:11px">No files attached</span>';
    } else { alert('Error'); }
}
async function updatePlasmid(e) {
    e.preventDefault();
    const fd = new FormData(e.target);
    const r = await fetch('/api/plasmids/'+fd.get('id'), {method:'PUT',body:fd});
    if (r.ok) location.reload(); else { const d = await r.json().catch(()=>({})); alert(d.error||'Error'); }
}

const featureColors = {'CDS':'#2e7d32','gene':'#2e7d32','promoter':'#e65100','terminator':'#c62828','rep_origin':'#1565c0','misc_feature':'#6a1b9a','protein_bind':'#00838f','RBS':'#9e9d24'};
async function viewSequence(filename) {
    document.getElementById('seqModal').classList.add('active');
    const dn = filename.includes('_')?filename.split('_').slice(3).join('_'):filename;
    document.getElementById('seqModalTitle').textContent = dn;
    document.getElementById('seqDownloadBtn').href = '/download/'+filename;
    document.getElementById('seqInfo').innerHTML = '<span style="color:var(--fg3)">Loading...</span>';
    document.getElementById('seqMapContainer').innerHTML = '<p style="color:var(--fg3);padding:40px">Loading...</p>';
    document.getElementById('seqText').textContent = 'Loading...';
    document.getElementById('seqFeaturesContainer').innerHTML = 'Loading...';
    showSeqTab('map');
    try {
        const resp = await fetch('/api/sequence/'+encodeURIComponent(filename));
        const data = await resp.json();
        let info = '<span class="seq-badge"><b>'+((data.length||0).toLocaleString())+' bp</b></span>';
        info += '<span class="seq-badge">'+( data.topology||'unknown')+'</span>';
        info += '<span class="seq-badge">'+( data.format||'unknown')+'</span>';
        if (data.features) info += '<span class="seq-badge">'+data.features.length+' features</span>';
        if (data.insert_region) info += '<span class="seq-badge highlight"><b>GOI:</b> '+data.insert_region.label+' ('+data.insert_region.start+'-'+data.insert_region.end+')</span>';
        document.getElementById('seqInfo').innerHTML = info;
        if (data.topology==='circular'&&data.length>0) {
            const sr = await fetch('/api/sequence-svg/'+encodeURIComponent(filename));
            document.getElementById('seqMapContainer').innerHTML = await sr.text();
        } else if (data.length>0) {
            document.getElementById('seqMapContainer').innerHTML = '<div style="padding:40px;color:var(--fg3)">Linear sequence ('+data.length.toLocaleString()+' bp)</div>';
        } else { document.getElementById('seqMapContainer').innerHTML = '<p style="color:var(--danger);padding:40px">Could not parse sequence data</p>'; }
        if (data.sequence&&!data.sequence.startsWith('[')) {
            let h=''; const s=data.sequence.toUpperCase(), ll=60, is=data.insert_region?data.insert_region.start:-1, ie=data.insert_region?data.insert_region.end:-1;
            for(let i=0;i<s.length;i+=ll){h+='<span class="line-num">'+String(i+1).padStart(6,' ')+'</span>';for(let j=i;j<Math.min(i+ll,s.length);j++){h+=j>=is&&j<ie?'<span class="insert-highlight">'+s[j]+'</span>':s[j];if((j-i+1)%10===0&&j<i+ll-1)h+=' ';}h+='\\n';}
            document.getElementById('seqText').innerHTML = h;
        } else { document.getElementById('seqText').textContent = data.sequence||'No sequence data'; }
        if (data.features&&data.features.length>0) {
            let ft='<table class="feature-table"><thead><tr><th>Type</th><th>Label</th><th>Location</th><th>Strand</th><th>Product</th></tr></thead><tbody>';
            data.features.forEach(f=>{const c=featureColors[f.type]||'#607D8B';const ins=data.insert_region&&f.start===data.insert_region.start&&f.end===data.insert_region.end;ft+='<tr'+(ins?' style="background:#fdf0ef;font-weight:600"':'')+'><td><span class="feature-type" style="background:'+c+'">'+f.type+'</span></td><td>'+(f.label||'-')+'</td><td>'+f.start+'..'+f.end+'</td><td>'+(f.strand===1?'+':f.strand===-1?'-':'?')+'</td><td>'+(f.product||'-')+'</td></tr>';});
            ft+='</tbody></table>';
            document.getElementById('seqFeaturesContainer').innerHTML = ft;
        } else { document.getElementById('seqFeaturesContainer').innerHTML = '<p style="color:var(--fg3);padding:20px">No features found.</p>'; }
    } catch(err) { document.getElementById('seqInfo').innerHTML = '<span style="color:var(--danger)">Error: '+err.message+'</span>'; }
}
function closeSeqModal() { document.getElementById('seqModal').classList.remove('active'); }
function showSeqTab(tab, evt) {
    document.querySelectorAll('.seq-tab-content').forEach(el=>el.style.display='none');
    document.querySelectorAll('.seq-tab').forEach(el=>el.classList.remove('active'));
    document.getElementById('seqTab'+tab.charAt(0).toUpperCase()+tab.slice(1)).style.display='block';
    if(evt&&evt.target){evt.target.classList.add('active');}else{const tabs=document.querySelectorAll('.seq-tab');const m={'map':0,'sequence':1,'features':2};if(tabs[m[tab]])tabs[m[tab]].classList.add('active');}
}
</script>
</body>
</html>
'''

@app.route('/')
def index():
    conn = sqlite3.connect(DB_FILE)
    c = conn.cursor()
    c.execute('SELECT * FROM plasmids ORDER BY number DESC')
    plasmids = c.fetchall()
    conn.close()
    return render_template_string(HTML_TEMPLATE, plasmids=plasmids)

@app.route('/uploads/<filename>')
def uploaded_file(filename):
    return send_from_directory(UPLOAD_DIR, filename)

@app.route('/download/<filename>')
def download_file(filename):
    return send_from_directory(UPLOAD_DIR, filename, as_attachment=True)

@app.route('/api/next-id', methods=['GET'])
def get_next_id():
    """Get the next available ID number"""
    conn = sqlite3.connect(DB_FILE)
    c = conn.cursor()
    c.execute('SELECT MAX(number) FROM plasmids')
    max_id = c.fetchone()[0]
    conn.close()
    next_id = (max_id or 0) + 1
    return jsonify({'next_id': next_id})

@app.route('/api/plasmids', methods=['POST'])
def add_plasmid():
    ip = request.remote_addr
    
    # Handle file uploads
    uploaded_files = []
    if 'files' in request.files:
        files = request.files.getlist('files')
        for file in files:
            if file and file.filename and allowed_file(file.filename):
                filename = secure_filename(file.filename)
                # Add timestamp to avoid conflicts
                timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
                filename = f"{timestamp}_{filename}"
                file.save(os.path.join(UPLOAD_DIR, filename))
                uploaded_files.append(filename)
    
    conn = sqlite3.connect(DB_FILE)
    c = conn.cursor()
    
    # Auto-increment ID if not provided
    number = request.form.get('number')
    if not number:
        c.execute('SELECT MAX(number) FROM plasmids')
        max_id = c.fetchone()[0]
        number = (max_id or 0) + 1
    
    files_str = ','.join(uploaded_files) if uploaded_files else None
    
    data = {
        'name': request.form.get('name'),
        'number': number,
        'plasmid_backbone': request.form.get('plasmid_backbone'),
        'inductor': request.form.get('inductor'),
        'antibiotic_resistance': request.form.get('antibiotic_resistance'),
        'reference': request.form.get('reference'),
        'comment': request.form.get('comment'),
        'files': files_str,
        'added_by': request.form.get('added_by')
    }
    
    # Validate mandatory fields
    if not data['added_by']:
        return jsonify({'error': 'Added By is required'}), 400
    c.execute('''INSERT INTO plasmids 
                 (name, number, plasmid_backbone, inductor, antibiotic_resistance, reference, comment, files, added_by)
                 VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)''',
              (data['name'], data['number'], data['plasmid_backbone'],
               data['inductor'], data['antibiotic_resistance'], 
               data['reference'], data['comment'], data['files'], data['added_by']))
    record_id = c.lastrowid
    conn.commit()
    conn.close()
    
    # Log to audit trail
    log_audit('INSERT', record_id=record_id, new_data=data, ip_address=ip)
    
    return jsonify({'success': True, 'id': record_id})

@app.route('/api/plasmids/<int:id>', methods=['PUT'])
def update_plasmid(id):
    ip = request.remote_addr
    
    conn = sqlite3.connect(DB_FILE)
    c = conn.cursor()
    
    # Get old data
    c.execute('SELECT * FROM plasmids WHERE id = ?', (id,))
    old_data = c.fetchone()
    
    # Handle file uploads
    uploaded_files = []
    if 'files' in request.files:
        files = request.files.getlist('files')
        for file in files:
            if file and file.filename and allowed_file(file.filename):
                filename = secure_filename(file.filename)
                timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
                filename = f"{timestamp}_{filename}"
                file.save(os.path.join(UPLOAD_DIR, filename))
                uploaded_files.append(filename)
    
    # Merge with existing files (column index 10 = files)
    existing_files = old_data[10] if old_data[10] else ''
    all_files = []
    if existing_files:
        all_files.extend(existing_files.split(','))
    all_files.extend(uploaded_files)
    files_str = ','.join(all_files) if all_files else None
    
    new_data = {
        'name': request.form.get('name'),
        'number': old_data[2],  # ID CANNOT CHANGE
        'plasmid_backbone': request.form.get('plasmid_backbone'),
        'inductor': request.form.get('inductor'),
        'antibiotic_resistance': request.form.get('antibiotic_resistance'),
        'reference': request.form.get('reference'),
        'comment': request.form.get('comment'),
        'files': files_str,
        'added_by': request.form.get('added_by') or old_data[12]
    }
    
    c.execute('''UPDATE plasmids 
                 SET name=?, plasmid_backbone=?, inductor=?, antibiotic_resistance=?, 
                     reference=?, comment=?, files=?, added_by=?, modified_at=CURRENT_TIMESTAMP
                 WHERE id=?''',
              (new_data['name'], new_data['plasmid_backbone'], new_data['inductor'],
               new_data['antibiotic_resistance'], new_data['reference'], 
               new_data['comment'], new_data['files'], new_data['added_by'], id))
    conn.commit()
    conn.close()
    
    # Log to audit trail
    log_audit('UPDATE', record_id=id, old_data=old_data, new_data=new_data, ip_address=ip)
    
    return jsonify({'success': True})

@app.route('/api/plasmids/<int:id>', methods=['DELETE'])
def delete_plasmid(id):
    """DELETE ONLY AVAILABLE IN ADMIN PANEL"""
    ip = request.remote_addr
    
    # Get old data before deleting
    conn = sqlite3.connect(DB_FILE)
    c = conn.cursor()
    c.execute('SELECT * FROM plasmids WHERE id = ?', (id,))
    old_data = c.fetchone()
    
    c.execute('DELETE FROM plasmids WHERE id = ?', (id,))
    conn.commit()
    conn.close()
    
    # Log to audit trail
    log_audit('DELETE', record_id=id, old_data=old_data, ip_address=ip)
    
    return jsonify({'success': True})

@app.route('/api/plasmids/<int:id>/remove-file', methods=['POST'])
def remove_file(id):
    """Remove a specific file from a plasmid entry"""
    ip = request.remote_addr
    filename = request.json.get('filename')
    if not filename:
        return jsonify({'error': 'No filename provided'}), 400

    conn = sqlite3.connect(DB_FILE)
    c = conn.cursor()
    c.execute('SELECT * FROM plasmids WHERE id = ?', (id,))
    old_data = c.fetchone()
    if not old_data:
        conn.close()
        return jsonify({'error': 'Entry not found'}), 404

    current_files = old_data[10] or ''
    file_list = [f.strip() for f in current_files.split(',') if f.strip()]
    if filename not in file_list:
        conn.close()
        return jsonify({'error': 'File not found on this entry'}), 404

    file_list.remove(filename)
    new_files = ','.join(file_list) if file_list else None

    c.execute('UPDATE plasmids SET files=?, modified_at=CURRENT_TIMESTAMP WHERE id=?',
              (new_files, id))
    conn.commit()
    conn.close()

    log_audit('REMOVE_FILE', record_id=id,
              old_data={'file': filename}, new_data={'files': new_files},
              ip_address=ip)

    return jsonify({'success': True, 'remaining_files': new_files})

@app.route('/api/backup', methods=['POST'])
def manual_backup():
    """Endpoint for manual backup"""
    backup_database()
    return jsonify({'success': True, 'message': 'Backup created'})

def parse_sequence_file(filepath):
    """Parse a sequence file and return sequence data + features"""
    ext = os.path.splitext(filepath)[1].lower()
    result = {
        'sequence': '',
        'features': [],
        'name': '',
        'length': 0,
        'topology': 'linear',
        'format': ext,
        'insert_region': None
    }
    
    try:
        if ext in ('.gb', '.gbk'):
            record = SeqIO.read(filepath, 'genbank')
            result['sequence'] = str(record.seq)
            result['name'] = record.description or record.name
            result['length'] = len(record.seq)
            if hasattr(record, 'annotations') and record.annotations.get('topology') == 'circular':
                result['topology'] = 'circular'
            elif 'circular' in open(filepath).readline().lower():
                result['topology'] = 'circular'
            for feat in record.features:
                f = {
                    'type': feat.type,
                    'start': int(feat.location.start),
                    'end': int(feat.location.end),
                    'strand': feat.location.strand,
                    'label': feat.qualifiers.get('label', feat.qualifiers.get('gene', feat.qualifiers.get('note', [''])))[0] if feat.qualifiers else '',
                    'product': feat.qualifiers.get('product', [''])[0] if feat.qualifiers else ''
                }
                result['features'].append(f)
                # Detect insert / gene of interest
                note = ' '.join(feat.qualifiers.get('note', []) + feat.qualifiers.get('label', []))
                if any(kw in note.lower() for kw in ['insert', 'gene of interest', 'goi']):
                    result['insert_region'] = {'start': f['start'], 'end': f['end'], 'label': f['label']}
                    
        elif ext == '.ape':
            # APE files are GenBank format
            record = SeqIO.read(filepath, 'genbank')
            result['sequence'] = str(record.seq)
            result['name'] = record.description or record.name
            result['length'] = len(record.seq)
            for feat in record.features:
                f = {
                    'type': feat.type,
                    'start': int(feat.location.start),
                    'end': int(feat.location.end),
                    'strand': feat.location.strand,
                    'label': feat.qualifiers.get('label', feat.qualifiers.get('gene', feat.qualifiers.get('note', [''])))[0] if feat.qualifiers else '',
                    'product': feat.qualifiers.get('product', [''])[0] if feat.qualifiers else ''
                }
                result['features'].append(f)
            # Check for circular topology
            with open(filepath) as fh:
                first_line = fh.readline()
                if 'circular' in first_line.lower():
                    result['topology'] = 'circular'
                    
        elif ext == '.dna':
            # SnapGene format
            try:
                record = SeqIO.read(filepath, 'snapgene')
                result['sequence'] = str(record.seq)
                result['name'] = record.name or os.path.basename(filepath)
                result['length'] = len(record.seq)
                result['topology'] = 'circular'  # Most SnapGene plasmids are circular
                for feat in record.features:
                    f = {
                        'type': feat.type,
                        'start': int(feat.location.start),
                        'end': int(feat.location.end),
                        'strand': feat.location.strand,
                        'label': feat.qualifiers.get('label', feat.qualifiers.get('gene', feat.qualifiers.get('note', [''])))[0] if feat.qualifiers else '',
                        'product': feat.qualifiers.get('product', [''])[0] if feat.qualifiers else ''
                    }
                    result['features'].append(f)
            except Exception:
                # Fallback: try reading raw bytes for SnapGene
                with open(filepath, 'rb') as fh:
                    data = fh.read()
                # SnapGene files start with byte 0x09
                if data[0:1] == b'\x09':
                    result['name'] = os.path.basename(filepath)
                    result['sequence'] = '[SnapGene binary - could not fully parse]'
                    
        elif ext == '.fasta':
            records = list(SeqIO.parse(filepath, 'fasta'))
            if records:
                record = records[0]
                result['sequence'] = str(record.seq)
                result['name'] = record.description or record.id
                result['length'] = len(record.seq)
                
        elif ext == '.seq':
            # .seq files can be FASTA or plain sequence
            with open(filepath) as fh:
                content = fh.read().strip()
            if content.startswith('>'):
                # FASTA format
                handle = io.StringIO(content)
                try:
                    record = next(SeqIO.parse(handle, 'fasta'))
                    result['sequence'] = str(record.seq)
                    result['name'] = record.description or record.id
                    result['length'] = len(record.seq)
                except StopIteration:
                    result['sequence'] = content
                    result['length'] = len(content)
            else:
                # Plain sequence text
                seq = re.sub(r'[^ATCGatcgNnRrYyMmKkSsWwBbDdHhVv]', '', content)
                result['sequence'] = seq
                result['name'] = os.path.basename(filepath)
                result['length'] = len(seq)
        else:
            result['sequence'] = '[Unsupported format]'
            
    except Exception as e:
        result['error'] = str(e)
        result['sequence'] = f'[Error parsing file: {e}]'
    
    return result

def generate_plasmid_svg(seq_data, width=500, height=500):
    """Generate a circular plasmid map as SVG"""
    cx, cy = width / 2, height / 2
    radius = min(width, height) * 0.35
    seq_len = seq_data['length'] or 1
    
    # Color palette for features
    colors = {
        'CDS': '#4CAF50', 'gene': '#4CAF50',
        'promoter': '#FF9800', 
        'terminator': '#f44336',
        'rep_origin': '#2196F3', 'ori': '#2196F3',
        'misc_feature': '#9C27B0',
        'protein_bind': '#00BCD4',
        'RBS': '#CDDC39',
        'default': '#607D8B'
    }
    
    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" style="max-width:100%;height:auto">
    <defs>
        <marker id="arrowhead" markerWidth="6" markerHeight="4" refX="3" refY="2" orient="auto">
            <polygon points="0 0, 6 2, 0 4" fill="#333"/>
        </marker>
    </defs>
    <style>
        text {{ font-family: -apple-system, sans-serif; }}
        .plasmid-name {{ font-size: 14px; font-weight: 600; fill: #2c3e50; text-anchor: middle; }}
        .plasmid-size {{ font-size: 11px; fill: #7f8c8d; text-anchor: middle; }}
        .feature-label {{ font-size: 8px; fill: #333; }}
        .tick-label {{ font-size: 7px; fill: #999; text-anchor: middle; }}
    </style>'''
    
    # Draw backbone circle
    svg += f'\n    <circle cx="{cx}" cy="{cy}" r="{radius}" fill="none" stroke="#bdc3c7" stroke-width="3"/>'
    
    # Draw tick marks every 1000bp
    for bp in range(0, seq_len, 1000):
        angle = (bp / seq_len) * 2 * math.pi - math.pi / 2
        x1 = cx + (radius - 8) * math.cos(angle)
        y1 = cy + (radius - 8) * math.sin(angle)
        x2 = cx + (radius + 8) * math.cos(angle)
        y2 = cy + (radius + 8) * math.sin(angle)
        xt = cx + (radius + 18) * math.cos(angle)
        yt = cy + (radius + 18) * math.sin(angle)
        svg += f'\n    <line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="#ccc" stroke-width="1"/>'
        if bp % 2000 == 0:
            svg += f'\n    <text x="{xt:.1f}" y="{yt:.1f}" class="tick-label">{bp//1000}k</text>'
    
    # Draw features as arcs
    insert_region = seq_data.get('insert_region')
    for feat in seq_data.get('features', []):
        start = feat['start']
        end = feat['end']
        label = feat['label'] or feat['type']
        ftype = feat['type']
        
        # Determine color
        color = colors.get(ftype, colors['default'])
        
        # Highlight insert region
        is_insert = False
        if insert_region and start == insert_region['start'] and end == insert_region['end']:
            is_insert = True
            color = '#e74c3c'
        
        # Calculate arc
        start_angle = (start / seq_len) * 2 * math.pi - math.pi / 2
        end_angle = (end / seq_len) * 2 * math.pi - math.pi / 2
        
        arc_r = radius + (4 if is_insert else 0)
        stroke_w = 10 if is_insert else 6
        
        x1 = cx + arc_r * math.cos(start_angle)
        y1 = cy + arc_r * math.sin(start_angle)
        x2 = cx + arc_r * math.cos(end_angle)
        y2 = cy + arc_r * math.sin(end_angle)
        
        large_arc = 1 if (end - start) / seq_len > 0.5 else 0
        
        svg += f'\n    <path d="M {x1:.1f} {y1:.1f} A {arc_r} {arc_r} 0 {large_arc} 1 {x2:.1f} {y2:.1f}" fill="none" stroke="{color}" stroke-width="{stroke_w}" opacity="0.8" stroke-linecap="round"/>'
        
        # Label
        mid_angle = (start_angle + end_angle) / 2
        if end_angle < start_angle:
            mid_angle += math.pi
        label_r = radius + 30 if is_insert else radius - 25
        lx = cx + label_r * math.cos(mid_angle)
        ly = cy + label_r * math.sin(mid_angle)
        
        # Only label features with names, skip very small ones
        if label and (end - start) > seq_len * 0.01:
            display_label = label[:20]
            fw = 'bold' if is_insert else 'normal'
            fs = '9px' if is_insert else '8px'
            svg += f'\n    <text x="{lx:.1f}" y="{ly:.1f}" class="feature-label" style="font-weight:{fw};font-size:{fs}" text-anchor="middle" dominant-baseline="middle">{display_label}</text>'
    
    # Center text
    name = (seq_data.get('name') or 'Plasmid')[:25]
    svg += f'\n    <text x="{cx}" y="{cy - 8}" class="plasmid-name">{name}</text>'
    svg += f'\n    <text x="{cx}" y="{cy + 12}" class="plasmid-size">{seq_len:,} bp</text>'
    
    svg += '\n</svg>'
    return svg

@app.route('/api/sequence/<filename>')
def get_sequence(filename):
    """Parse and return sequence data from a file"""
    filepath = os.path.join(UPLOAD_DIR, filename)
    if not os.path.exists(filepath):
        return jsonify({'error': 'File not found'}), 404
    
    seq_data = parse_sequence_file(filepath)
    return jsonify(seq_data)

@app.route('/api/sequence-svg/<filename>')
def get_sequence_svg(filename):
    """Return a circular plasmid map SVG"""
    filepath = os.path.join(UPLOAD_DIR, filename)
    if not os.path.exists(filepath):
        return Response('<svg><text>File not found</text></svg>', mimetype='image/svg+xml')
    
    seq_data = parse_sequence_file(filepath)
    if seq_data['length'] == 0:
        return Response('<svg><text>No sequence data</text></svg>', mimetype='image/svg+xml')
    
    svg = generate_plasmid_svg(seq_data)
    return Response(svg, mimetype='image/svg+xml')

@app.route('/export/html')
def export_html():
    """Generate a static HTML report of the entire plasmid library for Nagoya compliance"""
    conn = sqlite3.connect(DB_FILE)
    c = conn.cursor()
    c.execute('SELECT * FROM plasmids ORDER BY number ASC')
    plasmids = c.fetchall()
    conn.close()
    
    # Parse sequences for all plasmids that have files
    plasmid_data = []
    for p in plasmids:
        entry = {
            'id': p[0], 'name': p[1], 'number': p[2],
            'backbone': p[3], 'inductor': p[4],
            'resistance': p[5], 'reference': p[6],
            'comment': p[7], 'files': p[10],
            'created_at': p[9], 'added_by': p[12] if len(p) > 12 else '',
            'sequences': []
        }
        if p[10]:
            for fname in p[10].split(','):
                fname = fname.strip()
                filepath = os.path.join(UPLOAD_DIR, fname)
                ext = os.path.splitext(fname)[1].lower()
                if ext in ('.gb', '.gbk', '.ape', '.dna', '.fasta', '.seq') and os.path.exists(filepath):
                    seq_data = parse_sequence_file(filepath)
                    seq_data['filename'] = fname
                    if seq_data.get('sequence') and not seq_data['sequence'].startswith('['):
                        svg = generate_plasmid_svg(seq_data) if seq_data['topology'] == 'circular' and seq_data['length'] > 0 else ''
                        seq_data['svg'] = svg
                    entry['sequences'].append(seq_data)
        plasmid_data.append(entry)
    
    return render_template_string(EXPORT_TEMPLATE, 
                                  plasmids=plasmid_data,
                                  export_date=datetime.now().strftime('%Y-%m-%d %H:%M'),
                                  total=len(plasmids))

@app.route('/admin/download-backup/<filename>')
def download_backup(filename):
    """Download a specific backup file"""
    if not session.get('admin_logged_in'):
        return redirect('/admin/login')
    
    # Security: ensure filename is just a filename, not a path
    filename = os.path.basename(filename)
    backup_path = os.path.join(BACKUP_DIR, filename)
    
    if not os.path.exists(backup_path):
        return "Backup not found", 404
    
    return send_from_directory(BACKUP_DIR, filename, as_attachment=True)

@app.route('/admin/download-archive')
def download_archive():
    """Generate a ZIP archive containing the DB + self-contained HTML export"""
    if not session.get('admin_logged_in'):
        return redirect('/admin/login')
    
    # Generate the HTML export
    conn = sqlite3.connect(DB_FILE)
    c = conn.cursor()
    c.execute('SELECT * FROM plasmids ORDER BY number ASC')
    plasmids = c.fetchall()
    conn.close()
    
    plasmid_data = []
    for p in plasmids:
        entry = {
            'id': p[0], 'name': p[1], 'number': p[2],
            'backbone': p[3], 'inductor': p[4],
            'resistance': p[5], 'reference': p[6],
            'comment': p[7], 'files': p[10],
            'created_at': p[9], 'added_by': p[12] if len(p) > 12 else '',
            'sequences': []
        }
        if p[10]:
            for fname in p[10].split(','):
                fname = fname.strip()
                filepath = os.path.join(UPLOAD_DIR, fname)
                ext = os.path.splitext(fname)[1].lower()
                if ext in ('.gb', '.gbk', '.ape', '.dna', '.fasta', '.seq') and os.path.exists(filepath):
                    seq_data = parse_sequence_file(filepath)
                    seq_data['filename'] = fname
                    if seq_data.get('sequence') and not seq_data['sequence'].startswith('['):
                        svg = generate_plasmid_svg(seq_data) if seq_data['topology'] == 'circular' and seq_data['length'] > 0 else ''
                        seq_data['svg'] = svg
                    entry['sequences'].append(seq_data)
        plasmid_data.append(entry)
    
    export_date = datetime.now().strftime('%Y-%m-%d %H:%M')
    html_content = render_template_string(EXPORT_TEMPLATE, 
                                           plasmids=plasmid_data,
                                           export_date=export_date,
                                           total=len(plasmids))
    
    # Create ZIP in memory
    zip_buffer = io.BytesIO()
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zf:
        # Add database file
        zf.write(DB_FILE, f'plasmid_library_{timestamp}.db')
        # Add HTML export
        zf.writestr(f'plasmid_library_{timestamp}.html', html_content)
    
    zip_buffer.seek(0)
    return Response(
        zip_buffer.getvalue(),
        mimetype='application/zip',
        headers={'Content-Disposition': f'attachment; filename=plasmid_library_archive_{timestamp}.zip'}
    )

@app.route('/admin/login', methods=['GET', 'POST'])
def admin_login():
    """Admin login page"""
    if request.method == 'POST':
        password = request.form.get('password')
        if check_password_hash(ADMIN_PASSWORD_HASH, password):
            session['admin_logged_in'] = True
            return redirect('/admin')
        else:
            return render_template_string(LOGIN_TEMPLATE, error=True)
    return render_template_string(LOGIN_TEMPLATE, error=False)

@app.route('/admin/logout')
def admin_logout():
    """Admin logout"""
    session.pop('admin_logged_in', None)
    return redirect('/')

@app.route('/admin')
def admin_panel():
    """Admin panel with full database view and DELETE capability"""
    # Check if logged in
    if not session.get('admin_logged_in'):
        return redirect('/admin/login')
    
    conn = sqlite3.connect(DB_FILE)
    c = conn.cursor()
    
    # Get all plasmids
    c.execute('SELECT * FROM plasmids ORDER BY number DESC')
    plasmids = c.fetchall()
    
    # Get recent audit logs
    c.execute('SELECT * FROM audit_log ORDER BY id DESC LIMIT 100')
    audit_logs = c.fetchall()
    conn.close()
    
    # Get list of backups
    backups = []
    if os.path.exists(BACKUP_DIR):
        backup_files = sorted([f for f in os.listdir(BACKUP_DIR) if f.endswith('.db')], reverse=True)
        for f in backup_files:
            filepath = os.path.join(BACKUP_DIR, f)
            size = os.path.getsize(filepath) / 1024  # KB
            backups.append({
                'name': f,
                'size': f"{size:.1f} KB",
                'date': f.split('_')[3] + '_' + f.split('_')[4].replace('.db', '')
            })
    
    return render_template_string(ADMIN_TEMPLATE, plasmids=plasmids, audit_logs=audit_logs, backups=backups)

# Admin template
ADMIN_TEMPLATE = '''
<!DOCTYPE html>
<html><head><meta charset="UTF-8"><title>Admin - Plasmid Library</title>
<style>
:root { --fg: #1a1a1a; --fg2: #555; --fg3: #999; --bg: #fff; --bg2: #fafafa; --border: #e0e0e0; --accent: #1a6bb5; --danger: #c0392b; --mono: "SF Mono","Cascadia Code",Consolas,monospace; }
*,*::before,*::after{margin:0;padding:0;box-sizing:border-box}
html{font-size:13px}
body{font-family:"Helvetica Neue",Helvetica,Arial,sans-serif;color:var(--fg);background:var(--bg);line-height:1.5}
.topbar{display:flex;align-items:center;justify-content:space-between;padding:12px 24px;border-bottom:1px solid var(--border);position:sticky;top:0;z-index:100;background:var(--bg)}
.topbar h1{font-size:15px;font-weight:600}
.topbar-actions{display:flex;gap:6px}
.btn{display:inline-flex;align-items:center;gap:4px;padding:6px 14px;font-size:12px;font-weight:500;border:1px solid var(--border);border-radius:3px;background:var(--bg);color:var(--fg);cursor:pointer;text-decoration:none;white-space:nowrap}
.btn:hover{background:var(--bg2)}
.btn-primary{background:var(--fg);color:var(--bg);border-color:var(--fg)}
.btn-primary:hover{background:#333}
.btn-sm{padding:3px 8px;font-size:11px}
.btn-danger{color:var(--danger);border-color:var(--danger)}
.btn-danger:hover{background:#fdf0ef}
.content{padding:16px 24px}
.stats-row{display:flex;gap:24px;margin-bottom:20px;font-size:12px;color:var(--fg2)}
.stats-row strong{color:var(--fg);font-size:18px;display:block}
.section{margin-bottom:32px}
.section-head{display:flex;align-items:center;justify-content:space-between;padding:8px 0;border-bottom:2px solid var(--border);margin-bottom:12px;cursor:pointer;user-select:none}
.section-head h2{font-size:13px;font-weight:600;text-transform:uppercase;letter-spacing:0.03em;color:var(--fg2)}
.section-head .toggle{font-size:11px;color:var(--fg3)}
table{width:100%;border-collapse:collapse;font-size:12px}
th{padding:6px 8px;text-align:left;font-size:10px;font-weight:600;color:var(--fg2);text-transform:uppercase;letter-spacing:0.03em;border-bottom:2px solid var(--border);cursor:pointer;user-select:none;white-space:nowrap}
th:hover{color:var(--fg)}
td{padding:6px 8px;border-bottom:1px solid #f0f0f0;vertical-align:top}
tr:hover{background:#fafafa}
.badge{padding:2px 6px;font-size:10px;font-weight:600;font-family:var(--mono)}
.badge-insert{background:#e8f5e9;color:#2e7d32}
.badge-update{background:#fff8e1;color:#f57f17}
.badge-delete{background:#fce4ec;color:#c62828}
.badge-backup{background:#e3f2fd;color:#1565c0}
.backup-row{display:flex;justify-content:space-between;align-items:center;padding:8px 0;border-bottom:1px solid #f0f0f0;font-size:12px}
.backup-row:last-child{border-bottom:none}
/* modal */
.modal{display:none;position:fixed;z-index:1000;inset:0;background:rgba(0,0,0,0.4)}
.modal.active{display:flex;align-items:center;justify-content:center}
.modal-box{background:var(--bg);border:1px solid var(--border);width:90%;max-width:720px;max-height:90vh;overflow-y:auto}
.modal-head{display:flex;justify-content:space-between;align-items:center;padding:12px 16px;border-bottom:1px solid var(--border)}
.modal-head h2{font-size:14px;font-weight:600}
.modal-body{padding:16px}
.modal-foot{display:flex;gap:6px;padding:12px 16px;border-top:1px solid var(--border)}
.form-grid{display:grid;grid-template-columns:repeat(3,1fr);gap:12px}
.form-group{margin-bottom:8px}
.form-group.full-width{grid-column:1/-1}
.form-group label{display:block;font-size:11px;font-weight:600;color:var(--fg2);margin-bottom:4px;text-transform:uppercase;letter-spacing:0.03em}
.form-group input,.form-group textarea{width:100%;padding:6px 8px;font-size:12px;border:1px solid var(--border);border-radius:3px;font-family:inherit;background:var(--bg)}
.form-group input:focus,.form-group textarea:focus{outline:none;border-color:var(--accent)}
.form-group textarea{min-height:60px;resize:vertical}
.form-group input[readonly]{background:var(--bg2);color:var(--fg3)}
.file-tag{display:inline-flex;align-items:center;gap:6px;padding:4px 8px;margin:2px;background:var(--bg2);border:1px solid var(--border);font-size:11px}
.file-tag a{color:var(--accent);text-decoration:none}
.file-tag .remove-btn{background:none;border:none;color:var(--fg3);cursor:pointer;font-size:14px;line-height:1;padding:0 2px}
.file-tag .remove-btn:hover{color:var(--danger)}
</style></head>
<body>
<div class="topbar">
    <h1>Admin</h1>
    <div class="topbar-actions">
        <a href="/" class="btn">Back to Library</a>
        <a href="/admin/logout" class="btn">Logout</a>
    </div>
</div>
<div class="content">
<div class="stats-row">
    <div><strong>{{ plasmids|length }}</strong>Plasmids</div>
    <div><strong>{{ backups|length }}</strong>Backups</div>
</div>

<div class="section">
<div class="section-head" onclick="toggleSection('plasmids')"><h2>All Plasmids</h2><span class="toggle" id="plasmids-toggle">show</span></div>
<div id="plasmids-section" style="display:none">
<div style="overflow-x:auto"><table id="adminTable"><thead><tr>
<th onclick="sortAdminTable(0)">ID</th><th onclick="sortAdminTable(1)">Name</th><th onclick="sortAdminTable(2)">Backbone</th><th onclick="sortAdminTable(3)">Inductor</th><th onclick="sortAdminTable(4)">Resistance</th><th onclick="sortAdminTable(5)">Reference</th><th onclick="sortAdminTable(6)">Comment</th><th>Files</th><th onclick="sortAdminTable(8)">Added</th><th onclick="sortAdminTable(9)">By</th><th></th>
</tr></thead><tbody>
{% for p in plasmids %}<tr>
<td style="font-family:var(--mono);font-weight:600">{{ p[2] or '-' }}</td><td style="font-weight:500">{{ p[1] }}</td><td>{{ p[3] or '-' }}</td><td>{{ p[4] or '-' }}</td><td>{{ p[5] or '-' }}</td>
<td style="max-width:160px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap" title="{{ p[6] or '' }}">{{ p[6] or '-' }}</td>
<td style="max-width:200px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap" title="{{ p[7] or '' }}">{{ p[7] or '-' }}</td>
<td>{% if p[10] %}{{ p[10].split(',')|length }}{% else %}-{% endif %}</td>
<td>{{ p[11][:10] if p[11] else '-' }}</td><td>{{ p[12] or '-' }}</td>
<td style="white-space:nowrap"><button class="btn btn-sm" onclick='editPlasmid({{ p[0] }}, {{ p|tojson }})'>Edit</button> <button class="btn btn-sm btn-danger" onclick="deletePlasmid({{ p[0] }})">Delete</button></td>
</tr>{% endfor %}</tbody></table></div></div></div>

<div class="section">
<div class="section-head" onclick="toggleSection('backups')"><h2>Backups</h2><span class="toggle" id="backups-toggle">show</span></div>
<div id="backups-section" style="display:none">
<div style="margin-bottom:12px;display:flex;gap:6px"><button class="btn" onclick="createBackup()">Create Backup</button><a href="/admin/download-archive" class="btn">Download Archive</a></div>
{% if backups %}{% for b in backups %}<div class="backup-row"><div><span style="font-weight:500">{{ b.name }}</span><span style="color:var(--fg3);margin-left:8px">{{ b.size }}</span></div><a href="/admin/download-backup/{{ b.name }}" class="btn btn-sm">Download</a></div>{% endfor %}{% else %}<p style="color:var(--fg3)">No backups yet</p>{% endif %}
</div></div>

<div class="section">
<div class="section-head" onclick="toggleSection('audit')"><h2>Audit Log</h2><span class="toggle" id="audit-toggle">show</span></div>
<div id="audit-section" style="display:none">
{% if audit_logs %}<div style="overflow-x:auto"><table><thead><tr><th>ID</th><th>Time</th><th>Action</th><th>Record</th><th>IP</th><th>Data</th></tr></thead><tbody>
{% for log in audit_logs %}<tr><td>{{ log[0] }}</td><td style="font-family:var(--mono);font-size:11px">{{ log[1][:19] }}</td><td><span class="badge badge-{{ log[2].lower() }}">{{ log[2] }}</span></td><td>{{ log[3] or '-' }}</td><td style="font-family:var(--mono)">{{ log[6] or '-' }}</td><td style="max-width:240px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap;font-family:var(--mono);font-size:10px">{{ log[4][:40] if log[4] else (log[5][:40] if log[5] else '-') }}</td></tr>{% endfor %}
</tbody></table></div>{% else %}<p style="color:var(--fg3)">No logs</p>{% endif %}
</div></div>
</div>

<div id="editModal" class="modal">
<div class="modal-box">
<div class="modal-head"><h2>Edit Plasmid</h2><button class="btn btn-sm" onclick="closeEditModal()">Close</button></div>
<form onsubmit="updatePlasmid(event)" enctype="multipart/form-data">
<div class="modal-body">
<input type="hidden" id="edit_id" name="id">
<div class="form-grid">
<div class="form-group"><label>ID</label><input type="number" id="edit_number" name="number" readonly></div>
<div class="form-group"><label>Name *</label><input type="text" id="edit_name" name="name" required></div>
<div class="form-group"><label>Backbone</label><input type="text" id="edit_plasmid_backbone" name="plasmid_backbone"></div>
<div class="form-group"><label>Inductor</label><input type="text" id="edit_inductor" name="inductor"></div>
<div class="form-group"><label>Resistance</label><input type="text" id="edit_antibiotic_resistance" name="antibiotic_resistance"></div>
<div class="form-group"><label>Reference</label><input type="text" id="edit_reference" name="reference"></div>
<div class="form-group"><label>Added By</label><input type="text" id="edit_added_by" name="added_by"></div>
<div class="form-group full-width"><label>Comment</label><textarea id="edit_comment" name="comment"></textarea></div>
<div class="form-group full-width"><label>Files</label><div id="current_files"></div></div>
<div class="form-group full-width"><label>Add Files</label><input type="file" name="files" multiple></div>
</div>
</div>
<div class="modal-foot"><button type="submit" class="btn btn-primary">Save</button><button type="button" class="btn" onclick="closeEditModal()">Cancel</button></div>
</form>
</div>
</div>

<script>
let adminSortDir={}, _editId=null;
function toggleSection(id){const s=document.getElementById(id+'-section'),t=document.getElementById(id+'-toggle');if(s.style.display==='none'){s.style.display='block';t.textContent='hide'}else{s.style.display='none';t.textContent='show'}}
function sortAdminTable(col){const table=document.getElementById('adminTable'),tbody=table.querySelector('tbody'),rows=Array.from(tbody.querySelectorAll('tr'));const dir=adminSortDir[col]==='asc'?'desc':'asc';adminSortDir={[col]:dir};rows.sort((a,b)=>{let av=a.cells[col].textContent.trim(),bv=b.cells[col].textContent.trim();if(av==='-')av='';if(bv==='-')bv='';const an=parseFloat(av),bn=parseFloat(bv);if(!isNaN(an)&&!isNaN(bn))return dir==='asc'?an-bn:bn-an;return dir==='asc'?av.localeCompare(bv):bv.localeCompare(av)});rows.forEach(r=>tbody.appendChild(r))}
async function createBackup(){if(!confirm('Create backup?'))return;const r=await fetch('/api/backup',{method:'POST'});if(r.ok)location.reload();else alert('Error')}
async function deletePlasmid(id){if(!confirm('Permanently delete this entry?'))return;const r=await fetch('/api/plasmids/'+id,{method:'DELETE'});if(r.ok)location.reload();else alert('Error')}
function editPlasmid(id,data){
    _editId=id;
    document.getElementById('edit_id').value=id;
    document.getElementById('edit_number').value=data[2]||'';
    document.getElementById('edit_name').value=data[1]||'';
    document.getElementById('edit_plasmid_backbone').value=data[3]||'';
    document.getElementById('edit_inductor').value=data[4]||'';
    document.getElementById('edit_antibiotic_resistance').value=data[5]||'';
    document.getElementById('edit_reference').value=data[6]||'';
    document.getElementById('edit_comment').value=data[7]||'';
    document.getElementById('edit_added_by').value=data[12]||'';
    const fd=document.getElementById('current_files');
    if(data[10]){fd.innerHTML=data[10].split(',').map(f=>{const dn=f.includes('_')?f.split('_').slice(3).join('_'):f;const safe=f.replace(/[^a-zA-Z0-9]/g,'_');return'<span class="file-tag" id="afile-'+safe+'"><a href="/download/'+f+'" download>'+dn+'</a><button type="button" class="remove-btn" data-fn="'+f+'">&times;</button></span>'}).join('');fd.querySelectorAll('.remove-btn').forEach(b=>b.onclick=function(){removeFile(this.dataset.fn)})}
    else{fd.innerHTML='<span style="color:var(--fg3);font-size:11px">No files</span>'}
    document.getElementById('editModal').classList.add('active');
}
function closeEditModal(){document.getElementById('editModal').classList.remove('active')}
async function removeFile(fn){if(!confirm('Remove file?'))return;const r=await fetch('/api/plasmids/'+_editId+'/remove-file',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({filename:fn})});if(r.ok){const el=document.getElementById('afile-'+fn.replace(/[^a-zA-Z0-9]/g,'_'));if(el)el.remove();const fd=document.getElementById('current_files');if(!fd.querySelector('.file-tag'))fd.innerHTML='<span style="color:var(--fg3);font-size:11px">No files</span>'}else alert('Error')}
async function updatePlasmid(e){e.preventDefault();const fd=new FormData(e.target);const r=await fetch('/api/plasmids/'+fd.get('id'),{method:'PUT',body:fd});if(r.ok)location.reload();else{const d=await r.json().catch(()=>({}));alert(d.error||'Error')}}
</script></body></html>
'''

# Login template
LOGIN_TEMPLATE = '''
<!DOCTYPE html>
<html><head><meta charset="UTF-8"><title>Admin Login</title>
<style>
*{margin:0;padding:0;box-sizing:border-box}
body{font-family:"Helvetica Neue",Helvetica,Arial,sans-serif;background:#fff;min-height:100vh;display:flex;align-items:center;justify-content:center}
.login-box{width:100%;max-width:320px;padding:24px}
.login-box h1{font-size:15px;font-weight:600;margin-bottom:4px}
.login-box p{font-size:12px;color:#999;margin-bottom:20px}
.field{margin-bottom:16px}
.field label{display:block;font-size:11px;font-weight:600;color:#555;text-transform:uppercase;letter-spacing:0.03em;margin-bottom:4px}
.field input{width:100%;padding:8px 10px;font-size:13px;border:1px solid #e0e0e0;border-radius:3px;font-family:inherit}
.field input:focus{outline:none;border-color:#1a6bb5}
.submit-btn{width:100%;padding:8px;font-size:13px;font-weight:500;background:#1a1a1a;color:#fff;border:none;border-radius:3px;cursor:pointer}
.submit-btn:hover{background:#333}
.error{padding:8px;margin-bottom:12px;font-size:12px;color:#c0392b;border:1px solid #fce4ec;background:#fdf0ef}
.back{margin-top:16px;text-align:center}
.back a{font-size:12px;color:#999;text-decoration:none}
.back a:hover{color:#1a1a1a}
</style></head>
<body>
<div class="login-box">
<h1>Admin Login</h1>
<p>Plasmid Library</p>
{% if error %}<div class="error">Incorrect password.</div>{% endif %}
<form method="POST">
<div class="field"><label>Password</label><input type="password" name="password" required autofocus></div>
<button type="submit" class="submit-btn">Login</button>
</form>
<div class="back"><a href="/">Back to library</a></div>
</div>
</body></html>
'''

EXPORT_TEMPLATE = '''
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Plasmid Library Registry</title>
<style>
:root { --fg: #1a1a1a; --fg2: #555; --fg3: #999; --bg: #fff; --bg2: #fafafa; --border: #e0e0e0; --accent: #2e7d4f; --accent-light: #eaf6ef; --danger: #c0392b; --mono: "SF Mono", "Cascadia Code", "Fira Code", Consolas, monospace; }
*, *::before, *::after { margin: 0; padding: 0; box-sizing: border-box; }
body { font-family: "Helvetica Neue", Helvetica, Arial, sans-serif; color: var(--fg); background: var(--bg); line-height: 1.5; max-width: 1200px; margin: 0 auto; padding: 40px 30px; }

.hero { padding: 48px 0 40px; border-bottom: 3px solid var(--fg); margin-bottom: 36px; }
.hero h1 { font-size: 36px; font-weight: 700; letter-spacing: -0.02em; margin-bottom: 4px; color: var(--fg); }
.hero .subtitle { font-size: 15px; color: var(--fg2); margin-bottom: 20px; }
.hero-meta { display: flex; gap: 32px; flex-wrap: wrap; }
.hero-meta span { font-size: 13px; color: var(--fg2); }
.hero-meta strong { color: var(--fg); font-weight: 600; }

.toc { border: 1px solid var(--border); padding: 20px 24px; margin-bottom: 36px; }
.toc h2 { font-size: 13px; font-weight: 600; text-transform: uppercase; letter-spacing: 0.04em; color: var(--fg2); margin-bottom: 12px; }
.toc-list { column-count: 3; column-gap: 24px; list-style: none; font-size: 12px; }
.toc-list li { margin-bottom: 3px; }
.toc-list a { color: var(--accent); text-decoration: none; }
.toc-list a:hover { text-decoration: underline; }

.entry { margin-bottom: 28px; page-break-inside: avoid; border: 1px solid var(--border); border-left: 3px solid var(--accent); }
.entry-head { background: var(--fg); color: var(--bg); padding: 10px 20px; display: flex; justify-content: space-between; align-items: center; }
.entry-head h3 { font-size: 15px; font-weight: 600; }
.entry-head .pid { font-size: 22px; font-weight: 800; font-family: var(--mono); letter-spacing: -0.02em; color: #7dcea0; }
.entry-body { padding: 16px 20px; }
.entry-meta { display: grid; grid-template-columns: repeat(auto-fill, minmax(180px, 1fr)); gap: 8px; margin-bottom: 12px; }
.field { font-size: 12px; }
.field .label { display: block; font-size: 10px; font-weight: 600; text-transform: uppercase; letter-spacing: 0.03em; color: var(--fg3); }
.field .value { color: var(--fg); margin-top: 1px; }

.seq-section { margin-top: 14px; border-top: 1px solid var(--border); padding-top: 14px; }
.seq-section h4 { font-size: 13px; font-weight: 600; color: var(--fg); margin-bottom: 8px; }
.seq-badges { margin-bottom: 8px; }
.badge { display: inline-block; padding: 2px 8px; border: 1px solid var(--border); font-size: 11px; margin: 2px 2px 2px 0; background: var(--accent-light); }
.badge.insert { border-color: var(--danger); color: var(--danger); font-weight: 600; }
.seq-display { display: flex; gap: 20px; flex-wrap: wrap; }
.seq-map-export { flex: 0 0 280px; text-align: center; cursor: pointer; border: 1px solid var(--border); padding: 4px; }
.seq-map-export:hover { border-color: var(--fg); }
.seq-map-export::after { content: "click to zoom"; display: block; font-size: 10px; color: var(--fg3); margin-top: 4px; }
.seq-data-export { flex: 1; min-width: 280px; }
.seq-box { background: #1a1a1a; color: #ccc; font-family: var(--mono); font-size: 10px; padding: 12px; max-height: 280px; overflow-y: auto; word-break: break-all; white-space: pre-wrap; line-height: 1.6; }
.seq-box .insert-hl { color: #e74c3c; font-weight: bold; background: rgba(231,76,60,0.15); }

.features-mini { font-size: 11px; margin-top: 10px; }
.features-mini table { width: 100%; border-collapse: collapse; }
.features-mini th { padding: 4px 8px; text-align: left; font-size: 10px; font-weight: 600; text-transform: uppercase; letter-spacing: 0.03em; color: var(--fg2); border-bottom: 2px solid var(--fg); }
.features-mini td { padding: 4px 8px; border-bottom: 1px solid var(--border); font-size: 10px; }

.footer { border-top: 3px solid var(--accent); padding-top: 16px; margin-top: 48px; font-size: 12px; color: var(--fg3); text-align: center; }

.lightbox { display: none; position: fixed; z-index: 9999; inset: 0; background: rgba(0,0,0,0.85); align-items: center; justify-content: center; cursor: pointer; }
.lightbox.active { display: flex; }
.lightbox-content { max-width: 90vw; max-height: 90vh; background: white; padding: 20px; position: relative; }
.lightbox-content svg { width: 80vmin; height: 80vmin; }
.lightbox-close { position: absolute; top: -12px; right: -12px; background: var(--fg); color: white; border: none; width: 28px; height: 28px; font-size: 16px; cursor: pointer; font-weight: bold; }

@media print { body { padding: 10px; font-size: 11px; } .entry { break-inside: avoid; } .seq-box { max-height: none; } .lightbox { display: none !important; } }
</style>
</head>
<body>
<div class="hero">
<h1>Plasmid Library</h1>
<p class="subtitle">Genetic Material Registry</p>
<div class="hero-meta">
<span><strong>Export Date</strong> {{ export_date }}</span>
<span><strong>Total Entries</strong> {{ total }}</span>
<span><strong>Institution</strong> Your Institution</span>
</div>
</div>

<div class="toc">
<h2>Table of Contents ({{ total }} entries)</h2>
<ul class="toc-list">
{% for p in plasmids %}
<li><a href="#entry-{{ p.number }}">P{{ p.number }} {{ p.name }}</a></li>
{% endfor %}
</ul>
</div>

{% for p in plasmids %}
<div class="entry" id="entry-{{ p.number }}">
<div class="entry-head">
<h3>{{ p.name }}</h3>
<span class="pid">P{{ p.number }}</span>
</div>
<div class="entry-body">
<div class="entry-meta">
<div class="field"><span class="label">Backbone</span><span class="value">{{ p.backbone or '-' }}</span></div>
<div class="field"><span class="label">Inductor</span><span class="value">{{ p.inductor or '-' }}</span></div>
<div class="field"><span class="label">Resistance</span><span class="value">{{ p.resistance or '-' }}</span></div>
<div class="field"><span class="label">Reference</span><span class="value">{{ p.reference or '-' }}</span></div>
<div class="field"><span class="label">Date Added</span><span class="value">{{ p.created_at[:10] if p.created_at else '-' }}</span></div>
</div>
{% if p.comment %}<div class="field" style="margin-bottom:10px"><span class="label">Comment</span><span class="value">{{ p.comment }}</span></div>{% endif %}

{% for seq in p.sequences %}
<div class="seq-section">
<h4>{{ seq.filename.split('_', 3)[-1] if '_' in seq.filename else seq.filename }}</h4>
<div class="seq-badges">
<span class="badge"><strong>{{ seq.length|default(0) }}</strong> bp</span>
<span class="badge">{{ seq.topology|default('unknown') }}</span>
<span class="badge">{{ seq.format|default('unknown') }}</span>
{% if seq.insert_region %}<span class="badge insert"><strong>GOI:</strong> {{ seq.insert_region.label }} ({{ seq.insert_region.start }}-{{ seq.insert_region.end }})</span>{% endif %}
</div>
<div class="seq-display">
{% if seq.svg %}
<div class="seq-map-export" onclick="zoomMap(this)">{{ seq.svg|safe }}</div>
{% endif %}
<div class="seq-data-export">
{% if seq.sequence and not seq.sequence.startswith('[') %}
<div class="seq-box">{% if seq.insert_region %}{% set ins_s = seq.insert_region.start %}{% set ins_e = seq.insert_region.end %}{% for i in range(0, seq.sequence|length, 60) %}{% set line = seq.sequence[i:i+60] %}<span style="color:#666">{{ '%6d'|format(i+1) }} </span>{% for j in range(line|length) %}{% set pos = i + j %}{% if pos >= ins_s and pos < ins_e %}<span class="insert-hl">{{ line[j]|upper }}</span>{% else %}{{ line[j]|upper }}{% endif %}{% if (j+1) % 10 == 0 and j < line|length - 1 %} {% endif %}{% endfor %}
{% endfor %}{% else %}{% for i in range(0, seq.sequence|length, 60) %}{% set line = seq.sequence[i:i+60] %}<span style="color:#666">{{ '%6d'|format(i+1) }} </span>{% for j in range(line|length) %}{{ line[j]|upper }}{% if (j+1) % 10 == 0 and j < line|length - 1 %} {% endif %}{% endfor %}
{% endfor %}{% endif %}</div>
{% else %}
<p style="color:var(--fg3);font-size:12px">Sequence data could not be extracted from this file format.</p>
{% endif %}
</div>
</div>
{% if seq.features %}
<div class="features-mini">
<table><thead><tr><th>Type</th><th>Label</th><th>Location</th><th>Strand</th><th>Product</th></tr></thead><tbody>
{% for f in seq.features %}
<tr{% if seq.insert_region and f.start == seq.insert_region.start and f.end == seq.insert_region.end %} style="background:#fdf0ef;font-weight:600"{% endif %}>
<td>{{ f.type }}</td><td>{{ f.label or '-' }}</td><td>{{ f.start }}..{{ f.end }}</td><td>{{ '+' if f.strand == 1 else ('-' if f.strand == -1 else '?') }}</td><td>{{ f.product or '-' }}</td></tr>
{% endfor %}
</tbody></table>
</div>
{% endif %}
</div>
{% endfor %}
{% if not p.sequences %}<p style="color:var(--fg3);font-size:12px;margin-top:10px">No parseable sequence files attached.</p>{% endif %}
</div>
</div>
{% endfor %}

<div id="mapLightbox" class="lightbox" onclick="this.classList.remove('active')">
<div class="lightbox-content" onclick="event.stopPropagation()">
<button class="lightbox-close" onclick="document.getElementById('mapLightbox').classList.remove('active')">&times;</button>
<div id="lightboxSvg"></div>
</div>
</div>

<div class="footer">
<p>Plasmid Library Database | {{ export_date }} | Your Institution</p>
</div>

<script>
function zoomMap(el) {
    var svg = el.querySelector('svg');
    if (!svg) return;
    var clone = svg.cloneNode(true);
    clone.style.width = '80vmin';
    clone.style.height = '80vmin';
    document.getElementById('lightboxSvg').innerHTML = '';
    document.getElementById('lightboxSvg').appendChild(clone);
    document.getElementById('mapLightbox').classList.add('active');
}
</script>
</body>
</html>
'''

if __name__ == '__main__':
    init_db()
    
    # Start backup scheduler in background thread
    backup_thread = Thread(target=schedule_backups, daemon=True)
    backup_thread.start()
    
    print("=" * 60)
    print("🧬 Plasmid Library Server Starting...")
    print("=" * 60)
    print(f"Database file: {os.path.abspath(DB_FILE)}")
    print(f"Backup directory: {os.path.abspath(BACKUP_DIR)}")
    print(f"Upload directory: {os.path.abspath(UPLOAD_DIR)}")
    print(f"Main app: http://0.0.0.0:5001")
    print(f"Admin panel: http://0.0.0.0:5001/admin")
    print(f"Admin login: http://0.0.0.0:5001/admin/login")
    print(f"HTML Export: http://0.0.0.0:5001/export/html")
    print("=" * 60)
    print("🔐 ADMIN PASSWORD: admin123")
    print("⚠️  CHANGE THIS in plasmid_server.py line 33!")
    print("=" * 60)
    print("✓ Audit trail: Enabled (all changes tracked)")
    print("✓ Auto-backup: Every Friday at 6:00 PM")
    print("✓ File uploads: Enabled (50MB max)")
    print("✓ ID Protection: ID field is immutable")
    print("✓ Admin panel: Password protected")
    print("✓ Sequence viewer: Enabled (BioPython)")
    print("✓ Nagoya export: /export/html")
    print("=" * 60)
    
    app.run(host='0.0.0.0', port=5001, debug=False)
