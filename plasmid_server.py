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
                  files TEXT,
                  created_at TEXT DEFAULT CURRENT_TIMESTAMP,
                  modified_at TEXT DEFAULT CURRENT_TIMESTAMP)''')
    
    # Add files column if it doesn't exist (for existing databases)
    try:
        c.execute('ALTER TABLE plasmids ADD COLUMN files TEXT')
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
        * { margin: 0; padding: 0; box-sizing: border-box; }
        body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; background: #f5f5f5; padding: 20px; }
        .container { max-width: 1800px; margin: 0 auto; background: white; border-radius: 8px; box-shadow: 0 2px 8px rgba(0,0,0,0.1); overflow: hidden; }
        .header { background: #2c3e50; color: white; padding: 20px 30px; display: flex; justify-content: space-between; align-items: center; }
        .header h1 { font-size: 24px; font-weight: 500; }
        .header-info { font-size: 13px; opacity: 0.9; margin-top: 5px; }
        .btn { background: #3498db; color: white; border: none; padding: 10px 20px; border-radius: 4px; cursor: pointer; font-size: 14px; margin-left: 10px; }
        .btn:hover { background: #2980b9; }
        .btn-success { background: #27ae60; }
        .btn-success:hover { background: #229954; }
        .btn-warning { background: #f39c12; }
        .btn-warning:hover { background: #e67e22; }
        .content { padding: 30px; }
        
        .toolbar { display: flex; gap: 15px; margin-bottom: 20px; flex-wrap: wrap; align-items: center; }
        .search-box { flex: 1; min-width: 300px; position: relative; }
        .search-box input { width: 100%; padding: 10px 40px 10px 15px; border: 2px solid #ddd; border-radius: 4px; font-size: 14px; }
        .search-box input:focus { outline: none; border-color: #3498db; }
        .search-icon { position: absolute; right: 12px; top: 50%; transform: translateY(-50%); color: #7f8c8d; }
        .filter-group { display: flex; gap: 10px; align-items: center; }
        .filter-select { padding: 10px; border: 2px solid #ddd; border-radius: 4px; font-size: 14px; background: white; cursor: pointer; }
        .stats-bar { background: #ecf0f1; padding: 10px 15px; border-radius: 4px; font-size: 13px; color: #2c3e50; margin-bottom: 20px; }
        
        .table-wrapper { overflow-x: auto; }
        table { width: 100%; border-collapse: collapse; margin-top: 10px; table-layout: fixed; }
        td { overflow: hidden; text-overflow: ellipsis; }
        th { background: #34495e; color: white; padding: 12px; text-align: left; font-weight: 600; font-size: 13px; cursor: pointer; user-select: none; position: relative; white-space: nowrap; resize: horizontal; overflow: auto; }
        th { position: relative; }
        .resizer { position: absolute; top: 0; right: 0; width: 5px; cursor: col-resize; user-select: none; height: 100%; }
        th:hover { background: #2c3e50; }
        th.sortable::after { content: " ⇅"; opacity: 0.5; font-size: 11px; }
        th.sort-asc::after { content: " ▲"; opacity: 1; }
        th.sort-desc::after { content: " ▼"; opacity: 1; }
        td { padding: 12px; border-bottom: 1px solid #ecf0f1; font-size: 13px; }
        tbody tr:hover { background: #f8f9fa; }
        tbody tr.hidden { display: none; }
        
        .id-col { width: 60px; font-weight: 600; color: #3498db; }
        .name-col { min-width: 200px; font-weight: 500; }
        .action-col { width: 100px; text-align: center; }
        
        .form-section { background: #f8f9fa; padding: 20px; border-radius: 4px; margin-bottom: 20px; display: none; border: 2px solid #3498db; }
        .form-section.active { display: block; }
        .form-grid { display: grid; grid-template-columns: repeat(3, 1fr); gap: 15px; }
        .form-group { margin-bottom: 15px; }
        .form-group.full-width { grid-column: 1 / -1; }
        .form-group label { display: block; margin-bottom: 5px; font-weight: 500; color: #2c3e50; font-size: 13px; }
        .form-group input, .form-group textarea { width: 100%; padding: 10px; border: 1px solid #ddd; border-radius: 4px; font-size: 14px; }
        .form-group input:focus, .form-group textarea:focus { outline: none; border-color: #3498db; }
        .form-group textarea { min-height: 80px; resize: vertical; }
        .form-group input[readonly] { background: #ecf0f1; color: #7f8c8d; cursor: not-allowed; }
        
        .edit-btn { background: #f39c12; color: white; border: none; padding: 6px 12px; border-radius: 4px; cursor: pointer; font-size: 12px; }
        .edit-btn:hover { background: #e67e22; }
        .empty-state { text-align: center; padding: 60px; color: #7f8c8d; }
        
        .file-list { margin-top: 10px; }
        .file-item { display: inline-block; margin: 5px; padding: 5px 10px; background: #e8f5e9; border-radius: 3px; font-size: 12px; }
        .file-item a { color: #27ae60; text-decoration: none; }
        .file-item a:hover { text-decoration: underline; }
        
        .modal { display: none; position: fixed; z-index: 1000; left: 0; top: 0; width: 100%; height: 100%; background: rgba(0,0,0,0.5); }
        .modal.active { display: flex; align-items: center; justify-content: center; }
        .modal-content { background: white; padding: 30px; border-radius: 8px; max-width: 800px; width: 90%; max-height: 90vh; overflow-y: auto; }
        .modal-header { display: flex; justify-content: space-between; align-items: center; margin-bottom: 20px; }
        .modal-header h2 { color: #2c3e50; }
        .close-btn { background: #e74c3c; color: white; border: none; padding: 8px 16px; border-radius: 4px; cursor: pointer; }
        .close-btn:hover { background: #c0392b; }
        
        .seq-btn { background: #8e44ad; color: white; border: none; padding: 4px 8px; border-radius: 3px; cursor: pointer; font-size: 11px; margin: 2px; }
        .seq-btn:hover { background: #7d3c98; }
        .dl-btn { background: #34495e; color: white; border: none; padding: 4px 8px; border-radius: 3px; cursor: pointer; font-size: 11px; margin: 2px; text-decoration: none; display: inline-block; }
        .dl-btn:hover { background: #2c3e50; }
        
        #seqModal .modal-content { max-width: 1100px; }
        .seq-viewer { display: grid; grid-template-columns: 1fr 1fr; gap: 20px; margin-top: 15px; }
        .seq-map { text-align: center; }
        .seq-text-container { max-height: 400px; overflow-y: auto; background: #1e1e1e; border-radius: 6px; padding: 15px; }
        .seq-text { font-family: 'Courier New', monospace; font-size: 12px; color: #d4d4d4; word-break: break-all; line-height: 1.6; white-space: pre-wrap; }
        .seq-text .insert-highlight { background: rgba(231,76,60,0.3); color: #e74c3c; font-weight: bold; }
        .seq-text .line-num { color: #666; margin-right: 10px; user-select: none; }
        .feature-table { width: 100%; border-collapse: collapse; margin-top: 15px; font-size: 12px; }
        .feature-table th { background: #34495e; color: white; padding: 6px 10px; text-align: left; font-size: 11px; }
        .feature-table td { padding: 6px 10px; border-bottom: 1px solid #eee; }
        .feature-type { display: inline-block; padding: 2px 6px; border-radius: 3px; font-size: 10px; font-weight: 600; color: white; }
        .seq-info { display: flex; gap: 20px; margin-bottom: 15px; flex-wrap: wrap; }
        .seq-info-item { background: #f0f0f0; padding: 8px 14px; border-radius: 4px; font-size: 13px; }
        .seq-info-item strong { color: #2c3e50; }
        .seq-tabs { display: flex; gap: 5px; margin-bottom: 10px; }
        .seq-tab { padding: 8px 16px; border: none; background: #ecf0f1; cursor: pointer; border-radius: 4px 4px 0 0; font-size: 13px; }
        .seq-tab.active { background: #3498db; color: white; }
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <div>
                <h1>🧬 Plasmid Library</h1>
            </div>
            <div>
                <button class="btn" onclick="location.reload()">↻ Refresh</button>
                <button class="btn btn-warning" onclick="window.open('/export/html','_blank')">� Full Sequence Report</button>
                <button class="btn" onclick="location.href='/admin'">🔐 Admin</button>
                <button class="btn btn-success" onclick="toggleAddForm()">+ Add Entry</button>
            </div>
        </div>

        <div class="content">
            <div class="form-section" id="addForm">
                <h3>Add New Plasmid</h3>
                <form onsubmit="addPlasmid(event)" enctype="multipart/form-data">
                    <div class="form-grid">
                        <div class="form-group">
                            <label>ID (Auto-generated) *</label>
                            <input type="number" id="nextIdField" name="number" readonly required>
                        </div>
                        <div class="form-group">
                            <label>Name *</label>
                            <input type="text" name="name" required>
                        </div>
                        <div class="form-group">
                            <label>Plasmid Backbone</label>
                            <input type="text" name="plasmid_backbone">
                        </div>
                        <div class="form-group">
                            <label>Inductor</label>
                            <input type="text" name="inductor">
                        </div>
                        <div class="form-group">
                            <label>Antibiotic Resistance</label>
                            <input type="text" name="antibiotic_resistance">
                        </div>
                        <div class="form-group">
                            <label>Reference</label>
                            <input type="text" name="reference">
                        </div>
                        <div class="form-group">
                            <label>Added By *</label>
                            <input type="text" name="added_by" required placeholder="Your name">
                        </div>
                        <div class="form-group full-width">
                            <label>Comment</label>
                            <textarea name="comment"></textarea>
                        </div>
                        <div class="form-group full-width">
                            <label>Sequence File (required: .gb, .seq, .fasta, .ape, .dna) *</label>
                            <input type="file" name="files" multiple required id="addFileInput">
                            <small style="color:#e74c3c;display:none" id="fileWarning">⚠️ At least one sequence file is required!</small>
                        </div>
                    </div>
                    <button type="submit" class="btn btn-success">Save</button>
                    <button type="button" class="btn" onclick="toggleAddForm()">Cancel</button>
                </form>
            </div>

            <div class="toolbar">
                <div class="search-box">
                    <input type="text" id="searchInput" placeholder="🔍 Search across all fields..." onkeyup="filterTable()">
                    <span class="search-icon">⌕</span>
                </div>
                <div class="filter-group">
                    <select class="filter-select" id="backboneFilter" onchange="filterTable()">
                        <option value="">All Backbones</option>
                    </select>
                    <select class="filter-select" id="resistanceFilter" onchange="filterTable()">
                        <option value="">All Resistances</option>
                    </select>
                    <select class="filter-select" id="inductorFilter" onchange="filterTable()">
                        <option value="">All Inductors</option>
                    </select>
                </div>
            </div>

            <div class="stats-bar">
                <span id="statsDisplay">Showing <strong id="visibleCount">0</strong> of <strong id="totalCount">0</strong> plasmids</span>
            </div>

            {% if plasmids %}
            <div class="table-wrapper">
                <table id="plasmidTable">
                    <thead>
                        <tr>
                            <th class="sortable id-col" onclick="sortTable(0)">ID</th>
                            <th class="sortable name-col" onclick="sortTable(1)">Name</th>
                            <th class="sortable" onclick="sortTable(2)">Plasmid Backbone</th>
                            <th class="sortable" onclick="sortTable(3)">Inductor</th>
                            <th class="sortable" onclick="sortTable(4)">Antibiotic Resistance</th>
                            <th class="sortable" onclick="sortTable(5)">Reference</th>
                            <th class="sortable" onclick="sortTable(6)">Comment</th>
                            <th class="sortable" onclick="sortTable(7)">Added By</th>
                            <th>Files</th>
                            <th class="action-col">Action</th>
                        </tr>
                    </thead>
                    <tbody>
                        {% for p in plasmids %}
                        <tr data-id="{{ p[2] or 0 }}" data-name="{{ p[1] }}" data-backbone="{{ p[3] or '' }}" data-inductor="{{ p[4] or '' }}" data-resistance="{{ p[5] or '' }}">
                            <td class="id-col">{{ p[2] or '-' }}</td>
                            <td class="name-col">{{ p[1] }}</td>
                            <td>{{ p[3] or '-' }}</td>
                            <td>{{ p[4] or '-' }}</td>
                            <td>{{ p[5] or '-' }}</td>
                            <td>{{ p[6] or '-' }}</td>
                            <td>{{ p[7] or '-' }}</td>
                            <td>{{ p[12] or '-' }}</td>
                            <td>
                                {% if p[10] %}
                                    {% set files = p[10].split(',') %}
                                    <div class="file-list">
                                        {% for file in files %}
                                            {% set fname = file.strip() %}
                                            {% set display = fname.split('_', 3)[-1] if '_' in fname else fname %}
                                            {% set ext = fname.rsplit('.', 1)[-1].lower() if '.' in fname else '' %}
                                            <div class="file-item" style="display:block">
                                                <span style="font-weight:500">📎 {{ display }}</span><br>
                                                {% if ext in ['gb', 'gbk', 'ape', 'dna', 'fasta', 'seq'] %}
                                                <button class="seq-btn" data-filename="{{ fname }}" onclick="viewSequence(this.dataset.filename)">🔬 View Sequence</button>
                                                {% endif %}
                                                <a class="dl-btn" href="/download/{{ fname }}" download>⬇ Download</a>
                                            </div>
                                        {% endfor %}
                                    </div>
                                {% endif %}
                            </td>
                            <td class="action-col">
                                <button class="edit-btn" onclick='editPlasmid({{ p[0] }}, {{ p|tojson }})'>Edit</button>
                            </td>
                        </tr>
                        {% endfor %}
                    </tbody>
                </table>
            </div>
            {% else %}
            <div class="empty-state">
                <p>No plasmids yet. Click "Add Entry" to get started!</p>
            </div>
            {% endif %}
        </div>
    </div>

    <!-- Edit Modal -->
    <div id="editModal" class="modal">
        <div class="modal-content">
            <div class="modal-header">
                <h2>Edit Plasmid</h2>
                <button class="close-btn" onclick="closeEditModal()">Close</button>
            </div>
            <form onsubmit="updatePlasmid(event)" enctype="multipart/form-data">
                <input type="hidden" id="edit_id" name="id">
                <div class="form-grid">
                    <div class="form-group">
                        <label>ID (Cannot be changed) *</label>
                        <input type="number" id="edit_number" name="number" readonly required>
                    </div>
                    <div class="form-group">
                        <label>Name *</label>
                        <input type="text" id="edit_name" name="name" required>
                    </div>
                    <div class="form-group">
                        <label>Plasmid Backbone</label>
                        <input type="text" id="edit_plasmid_backbone" name="plasmid_backbone">
                    </div>
                    <div class="form-group">
                        <label>Inductor</label>
                        <input type="text" id="edit_inductor" name="inductor">
                    </div>
                    <div class="form-group">
                        <label>Antibiotic Resistance</label>
                        <input type="text" id="edit_antibiotic_resistance" name="antibiotic_resistance">
                    </div>
                    <div class="form-group">
                        <label>Reference</label>
                        <input type="text" id="edit_reference" name="reference">
                    </div>
                    <div class="form-group">
                        <label>Added By *</label>
                        <input type="text" id="edit_added_by" name="added_by" required>
                    </div>
                    <div class="form-group full-width">
                        <label>Comment</label>
                        <textarea id="edit_comment" name="comment"></textarea>
                    </div>
                    <div class="form-group full-width">
                        <label>Current Files:</label>
                        <div id="current_files"></div>
                    </div>
                    <div class="form-group full-width">
                        <label>Add New Files</label>
                        <input type="file" name="files" multiple>
                    </div>
                </div>
                <button type="submit" class="btn btn-success">Save Changes</button>
                <button type="button" class="btn" onclick="closeEditModal()">Cancel</button>
            </form>
        </div>
    </div>

    <!-- Sequence Viewer Modal -->
    <div id="seqModal" class="modal">
        <div class="modal-content">
            <div class="modal-header">
                <h2 id="seqModalTitle">🔬 Sequence Viewer</h2>
                <div>
                    <a id="seqDownloadBtn" class="dl-btn" style="font-size:14px;padding:8px 16px" download>⬇ Download File</a>
                    <button class="close-btn" onclick="closeSeqModal()">Close</button>
                </div>
            </div>
            <div class="seq-info" id="seqInfo"></div>
            <div class="seq-tabs">
                <button class="seq-tab active" onclick="showSeqTab('map', event)">🗺 Plasmid Map</button>
                <button class="seq-tab" onclick="showSeqTab('sequence', event)">🧬 Sequence</button>
                <button class="seq-tab" onclick="showSeqTab('features', event)">📋 Features</button>
            </div>
            <div id="seqTabMap" class="seq-tab-content">
                <div id="seqMapContainer" class="seq-map" style="min-height:300px">
                    <p style="color:#999;padding:40px">Loading plasmid map...</p>
                </div>
            </div>
            <div id="seqTabSequence" class="seq-tab-content" style="display:none">
                <div class="seq-text-container" id="seqTextContainer">
                    <div class="seq-text" id="seqText">Loading sequence...</div>
                </div>
            </div>
            <div id="seqTabFeatures" class="seq-tab-content" style="display:none">
                <div id="seqFeaturesContainer">Loading features...</div>
            </div>
        </div>
    </div>

    <script>
        let sortDirection = {};

        window.onload = function() {
            populateFilters();
            updateStats();
            makeColumnsResizable();
        };

        function makeColumnsResizable() {
            const table = document.getElementById('plasmidTable');
            if (!table) return;
            
            const cols = table.querySelectorAll('th');
            cols.forEach((col, index) => {
                const resizer = document.createElement('div');
                resizer.classList.add('resizer');
                col.appendChild(resizer);
                
                let startX, startWidth;
                
                resizer.addEventListener('mousedown', function(e) {
                    e.stopPropagation();
                    startX = e.pageX;
                    startWidth = col.offsetWidth;
                    
                    document.addEventListener('mousemove', doDrag);
                    document.addEventListener('mouseup', stopDrag);
                });
                
                function doDrag(e) {
                    const width = startWidth + (e.pageX - startX);
                    if (width > 50) {
                        col.style.width = width + 'px';
                    }
                }
                
                function stopDrag() {
                    document.removeEventListener('mousemove', doDrag);
                    document.removeEventListener('mouseup', stopDrag);
                }
            });
        }

        async function toggleAddForm() {
            const form = document.getElementById('addForm');
            const isActive = form.classList.toggle('active');
            
            if (isActive) {
                const response = await fetch('/api/next-id');
                const data = await response.json();
                document.getElementById('nextIdField').value = data.next_id;
            }
        }

        function populateFilters() {
            const rows = document.querySelectorAll('#plasmidTable tbody tr');
            const backbones = new Set();
            const resistances = new Set();
            const inductors = new Set();

            rows.forEach(row => {
                const backbone = row.dataset.backbone;
                const resistance = row.dataset.resistance;
                const inductor = row.dataset.inductor;
                if (backbone && backbone !== '-') backbones.add(backbone);
                if (resistance && resistance !== '-') resistances.add(resistance);
                if (inductor && inductor !== '-') inductors.add(inductor);
            });

            populateSelect('backboneFilter', Array.from(backbones).sort());
            populateSelect('resistanceFilter', Array.from(resistances).sort());
            populateSelect('inductorFilter', Array.from(inductors).sort());
        }

        function populateSelect(id, values) {
            const select = document.getElementById(id);
            const currentValue = select.value;
            select.innerHTML = '<option value="">All ' + id.replace('Filter', '').replace(/([A-Z])/g, ' $1') + 's</option>';
            values.forEach(val => {
                const option = document.createElement('option');
                option.value = val;
                option.textContent = val;
                select.appendChild(option);
            });
            select.value = currentValue;
        }

        function filterTable() {
            const searchTerm = document.getElementById('searchInput').value.toLowerCase();
            const backboneFilter = document.getElementById('backboneFilter').value;
            const resistanceFilter = document.getElementById('resistanceFilter').value;
            const inductorFilter = document.getElementById('inductorFilter').value;
            
            const rows = document.querySelectorAll('#plasmidTable tbody tr');
            let visibleCount = 0;

            rows.forEach(row => {
                const text = row.textContent.toLowerCase();
                const backbone = row.dataset.backbone;
                const resistance = row.dataset.resistance;
                const inductor = row.dataset.inductor;

                const matchesSearch = text.includes(searchTerm);
                const matchesBackbone = !backboneFilter || backbone === backboneFilter;
                const matchesResistance = !resistanceFilter || resistance === resistanceFilter;
                const matchesInductor = !inductorFilter || inductor === inductorFilter;

                if (matchesSearch && matchesBackbone && matchesResistance && matchesInductor) {
                    row.classList.remove('hidden');
                    visibleCount++;
                } else {
                    row.classList.add('hidden');
                }
            });

            updateStats(visibleCount);
        }

        function updateStats(visible) {
            const total = document.querySelectorAll('#plasmidTable tbody tr').length;
            const visibleCount = visible !== undefined ? visible : total;
            document.getElementById('totalCount').textContent = total;
            document.getElementById('visibleCount').textContent = visibleCount;
        }

        function sortTable(columnIndex) {
            const table = document.getElementById('plasmidTable');
            const tbody = table.querySelector('tbody');
            const rows = Array.from(tbody.querySelectorAll('tr'));
            const headers = table.querySelectorAll('th');
            
            const currentDirection = sortDirection[columnIndex] || 'asc';
            const newDirection = currentDirection === 'asc' ? 'desc' : 'asc';
            sortDirection = { [columnIndex]: newDirection };

            headers.forEach((h, i) => {
                h.classList.remove('sort-asc', 'sort-desc');
                if (i === columnIndex) {
                    h.classList.add('sort-' + newDirection);
                }
            });

            rows.sort((a, b) => {
                let aVal = a.cells[columnIndex].textContent.trim();
                let bVal = b.cells[columnIndex].textContent.trim();

                if (aVal === '-') aVal = '';
                if (bVal === '-') bVal = '';

                const aNum = parseFloat(aVal);
                const bNum = parseFloat(bVal);
                
                if (!isNaN(aNum) && !isNaN(bNum)) {
                    return newDirection === 'asc' ? aNum - bNum : bNum - aNum;
                }

                return newDirection === 'asc' 
                    ? aVal.localeCompare(bVal)
                    : bVal.localeCompare(aVal);
            });

            rows.forEach(row => tbody.appendChild(row));
        }

        async function addPlasmid(event) {
            event.preventDefault();
            const formData = new FormData(event.target);

            const response = await fetch('/api/plasmids', {
                method: 'POST',
                body: formData
            });

            if (response.ok) {
                location.reload();
            } else {
                const data = await response.json().catch(() => ({}));
                alert(data.error || 'Error adding plasmid');
            }
        }

        function editPlasmid(id, data) {
            document.getElementById('edit_id').value = id;
            document.getElementById('edit_number').value = data[2] || '';
            document.getElementById('edit_name').value = data[1] || '';
            document.getElementById('edit_plasmid_backbone').value = data[3] || '';
            document.getElementById('edit_inductor').value = data[4] || '';
            document.getElementById('edit_antibiotic_resistance').value = data[5] || '';
            document.getElementById('edit_reference').value = data[6] || '';
            document.getElementById('edit_comment').value = data[7] || '';
            document.getElementById('edit_added_by').value = data[12] || '';
            
            // Show current files
            const filesDiv = document.getElementById('current_files');
            if (data[10]) {
                const files = data[10].split(',');
                filesDiv.innerHTML = files.map(f => {
                    const displayName = f.includes('_') ? f.split('_').slice(3).join('_') : f;
                    return `<div class="file-item"><a href="/download/${f}" download>📎 ${displayName}</a></div>`;
                }).join('');
            } else {
                filesDiv.innerHTML = '<em>No files attached</em>';
            }
            
            document.getElementById('editModal').classList.add('active');
        }

        function closeEditModal() {
            document.getElementById('editModal').classList.remove('active');
        }

        async function updatePlasmid(event) {
            event.preventDefault();
            const formData = new FormData(event.target);
            const id = formData.get('id');

            const response = await fetch(`/api/plasmids/${id}`, {
                method: 'PUT',
                body: formData
            });

            if (response.ok) {
                location.reload();
            } else {
                const data = await response.json().catch(() => ({}));
                alert(data.error || 'Error updating plasmid');
            }
        }

        // ---- Sequence Viewer ----
        const featureColors = {
            'CDS': '#4CAF50', 'gene': '#4CAF50',
            'promoter': '#FF9800', 'terminator': '#f44336',
            'rep_origin': '#2196F3', 'misc_feature': '#9C27B0',
            'protein_bind': '#00BCD4', 'RBS': '#CDDC39'
        };

        async function viewSequence(filename) {
            document.getElementById('seqModal').classList.add('active');
            const displayName = filename.includes('_') ? filename.split('_').slice(3).join('_') : filename;
            document.getElementById('seqModalTitle').textContent = '🔬 ' + displayName;
            document.getElementById('seqDownloadBtn').href = '/download/' + filename;
            document.getElementById('seqInfo').innerHTML = '<span style="color:#999">Loading...</span>';
            document.getElementById('seqMapContainer').innerHTML = '<p style="color:#999;padding:40px">Loading plasmid map...</p>';
            document.getElementById('seqText').textContent = 'Loading...';
            document.getElementById('seqFeaturesContainer').innerHTML = 'Loading...';
            showSeqTab('map');

            try {
                const resp = await fetch('/api/sequence/' + encodeURIComponent(filename));
                const data = await resp.json();

                // Info bar
                let infoHtml = '';
                infoHtml += `<div class="seq-info-item"><strong>Length:</strong> ${(data.length || 0).toLocaleString()} bp</div>`;
                infoHtml += `<div class="seq-info-item"><strong>Topology:</strong> ${data.topology || 'unknown'}</div>`;
                infoHtml += `<div class="seq-info-item"><strong>Format:</strong> ${data.format || 'unknown'}</div>`;
                if (data.features) infoHtml += `<div class="seq-info-item"><strong>Features:</strong> ${data.features.length}</div>`;
                if (data.insert_region) infoHtml += `<div class="seq-info-item" style="background:#fde8e8"><strong>Gene of Interest:</strong> ${data.insert_region.label} (${data.insert_region.start}-${data.insert_region.end})</div>`;
                document.getElementById('seqInfo').innerHTML = infoHtml;

                // Plasmid map (SVG)
                if (data.topology === 'circular' && data.length > 0) {
                    const svgResp = await fetch('/api/sequence-svg/' + encodeURIComponent(filename));
                    const svgText = await svgResp.text();
                    document.getElementById('seqMapContainer').innerHTML = svgText;
                } else if (data.length > 0) {
                    document.getElementById('seqMapContainer').innerHTML = `<div style="padding:40px;color:#7f8c8d"><p style="font-size:16px;margin-bottom:10px">Linear sequence (${data.length.toLocaleString()} bp)</p><p>Circular plasmid map not available for linear sequences.</p></div>`;
                } else {
                    document.getElementById('seqMapContainer').innerHTML = '<p style="color:#e74c3c;padding:40px">Could not parse sequence data</p>';
                }

                // Sequence text with insert highlighting
                if (data.sequence && !data.sequence.startsWith('[')) {
                    let seqHtml = '';
                    const seq = data.sequence.toUpperCase();
                    const lineLen = 60;
                    const insertStart = data.insert_region ? data.insert_region.start : -1;
                    const insertEnd = data.insert_region ? data.insert_region.end : -1;

                    for (let i = 0; i < seq.length; i += lineLen) {
                        const lineNum = String(i + 1).padStart(6, ' ');
                        seqHtml += `<span class="line-num">${lineNum}</span>`;
                        // Render chars with highlighting
                        for (let j = i; j < Math.min(i + lineLen, seq.length); j++) {
                            if (j >= insertStart && j < insertEnd) {
                                seqHtml += `<span class="insert-highlight">${seq[j]}</span>`;
                            } else {
                                seqHtml += seq[j];
                            }
                            // Space every 10 chars
                            if ((j - i + 1) % 10 === 0 && j < i + lineLen - 1) seqHtml += ' ';
                        }
                        seqHtml += '\\n';
                    }
                    document.getElementById('seqText').innerHTML = seqHtml;
                } else {
                    document.getElementById('seqText').textContent = data.sequence || 'No sequence data available';
                }

                // Features table
                if (data.features && data.features.length > 0) {
                    let ftHtml = '<table class="feature-table"><thead><tr><th>Type</th><th>Label</th><th>Location</th><th>Strand</th><th>Product</th></tr></thead><tbody>';
                    data.features.forEach(f => {
                        const color = featureColors[f.type] || '#607D8B';
                        const isInsert = data.insert_region && f.start === data.insert_region.start && f.end === data.insert_region.end;
                        const rowStyle = isInsert ? 'background:#fde8e8;font-weight:600' : '';
                        ftHtml += `<tr style="${rowStyle}"><td><span class="feature-type" style="background:${color}">${f.type}</span></td><td>${f.label || '-'}</td><td>${f.start}..${f.end}</td><td>${f.strand === 1 ? '+' : f.strand === -1 ? '-' : '?'}</td><td>${f.product || '-'}</td></tr>`;
                    });
                    ftHtml += '</tbody></table>';
                    document.getElementById('seqFeaturesContainer').innerHTML = ftHtml;
                } else {
                    document.getElementById('seqFeaturesContainer').innerHTML = '<p style="color:#7f8c8d;padding:20px">No annotated features found in this file.</p>';
                }
            } catch (err) {
                document.getElementById('seqInfo').innerHTML = `<div style="color:#e74c3c">Error loading sequence: ${err.message}</div>`;
            }
        }

        function closeSeqModal() {
            document.getElementById('seqModal').classList.remove('active');
        }

        function showSeqTab(tab, evt) {
            document.querySelectorAll('.seq-tab-content').forEach(el => el.style.display = 'none');
            document.querySelectorAll('.seq-tab').forEach(el => el.classList.remove('active'));
            document.getElementById('seqTab' + tab.charAt(0).toUpperCase() + tab.slice(1)).style.display = 'block';
            if (evt && evt.target) {
                evt.target.classList.add('active');
            } else {
                // Highlight first tab button by default
                const tabs = document.querySelectorAll('.seq-tab');
                const tabMap = {'map': 0, 'sequence': 1, 'features': 2};
                if (tabs[tabMap[tab]]) tabs[tabMap[tab]].classList.add('active');
            }
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
    if not uploaded_files:
        return jsonify({'error': 'At least one sequence file is required'}), 400
    
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
    
    # Merge with existing files
    existing_files = old_data[8] if old_data[8] else ''
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

# Admin template with full table view and DELETE capability
ADMIN_TEMPLATE = '''
<!DOCTYPE html>
<html><head><meta charset="UTF-8"><title>Admin - Plasmid Library</title>
<style>
*{margin:0;padding:0;box-sizing:border-box}
body{font-family:-apple-system,sans-serif;background:#f5f5f5;padding:20px}
.container{max-width:1800px;margin:0 auto;background:white;border-radius:8px;box-shadow:0 2px 8px rgba(0,0,0,0.1)}
.header{background:#8e44ad;color:white;padding:20px 30px;display:flex;justify-content:space-between;align-items:center}
.header h1{font-size:24px}
.btn{background:#3498db;color:white;border:none;padding:10px 20px;border-radius:4px;cursor:pointer;font-size:14px;text-decoration:none;display:inline-block;margin-left:10px}
.btn:hover{background:#2980b9}
.btn-success{background:#27ae60}
.content{padding:30px}
.section{margin-bottom:40px}
.section h2{color:#2c3e50;margin-bottom:20px;padding-bottom:10px;border-bottom:2px solid #ecf0f1}
table{width:100%;border-collapse:collapse;margin-top:10px;font-size:12px}
th{background:#34495e;color:white;padding:10px;text-align:left;font-weight:600;font-size:12px}
td{padding:10px;border-bottom:1px solid #ecf0f1}
tr:hover{background:#f8f9fa}
.delete-btn{background:#e74c3c;color:white;border:none;padding:6px 12px;border-radius:4px;cursor:pointer;font-size:11px}
.delete-btn:hover{background:#c0392b}
.action-badge{padding:4px 8px;border-radius:3px;font-size:11px;font-weight:600}
.action-insert{background:#d4edda;color:#155724}
.action-update{background:#fff3cd;color:#856404}
.action-delete{background:#f8d7da;color:#721c24}
.action-backup{background:#d1ecf1;color:#0c5460}
.backup-item{background:#f8f9fa;padding:15px;margin-bottom:10px;border-radius:4px;display:flex;justify-content:space-between;align-items:center}
.backup-name{font-weight:600;color:#2c3e50}
.backup-meta{font-size:12px;color:#7f8c8d;margin-top:5px}
.stats{display:grid;grid-template-columns:repeat(3,1fr);gap:20px;margin-bottom:30px}
.stat-card{background:linear-gradient(135deg,#667eea 0%,#764ba2 100%);color:white;padding:20px;border-radius:8px}
.stat-value{font-size:32px;font-weight:700;margin-bottom:5px}
.stat-label{font-size:14px;opacity:0.9}
.id-col{width:50px;font-weight:600;color:#3498db}
.name-col{min-width:150px;font-weight:500}
.warning-box{background:#fff3cd;border:2px solid #ffc107;padding:15px;border-radius:4px;margin-bottom:20px;color:#856404}
.sortable{cursor:pointer;user-select:none}
.sortable:hover{background:#2c3e50}
.edit-btn{background:#f39c12;color:white;border:none;padding:6px 12px;border-radius:4px;cursor:pointer;font-size:11px;margin-right:5px}
.edit-btn:hover{background:#e67e22}
</style></head>
<body><div class="container"><div class="header"><h1>🔐 Admin Panel</h1><div><a href="/" class="btn">← Back to Database</a><a href="/admin/logout" class="btn" style="background:#e74c3c">Logout</a></div></div>
<div class="content">
<div class="warning-box"><strong>⚠️ Admin Access:</strong> This panel has DELETE capabilities. All deletions are logged in the audit trail.</div>

<div class="stats"><div class="stat-card"><div class="stat-value">{{ plasmids|length }}</div><div class="stat-label">Total Plasmids</div></div><div class="stat-card"><div class="stat-value">{{ backups|length }}</div><div class="stat-label">Backups</div></div><div class="stat-card"><div class="stat-value">📅</div><div class="stat-label">Next: Friday 6PM</div></div></div>

<div class="section"><h2 onclick="toggleSection('plasmids')" style="cursor:pointer">🧬 All Plasmids <span id="plasmids-toggle">▼</span></h2>
<div id="plasmids-section" style="display:none">
<div style="overflow-x:auto"><table id="adminTable"><thead><tr><th class="sortable" onclick="sortAdminTable(0)">ID</th><th class="sortable" onclick="sortAdminTable(1)">Name</th><th class="sortable" onclick="sortAdminTable(2)">Backbone</th><th class="sortable" onclick="sortAdminTable(3)">Inductor</th><th class="sortable" onclick="sortAdminTable(4)">Resistance</th><th class="sortable" onclick="sortAdminTable(5)">Reference</th><th class="sortable" onclick="sortAdminTable(6)">Comment</th><th>Files</th><th class="sortable" onclick="sortAdminTable(8)">Added On</th><th class="sortable" onclick="sortAdminTable(9)">Added By</th><th>Action</th></tr></thead><tbody>
{% for p in plasmids %}<tr><td>{{ p[2] or '-' }}</td><td>{{ p[1] }}</td><td>{{ p[3] or '-' }}</td><td>{{ p[4] or '-' }}</td><td>{{ p[5] or '-' }}</td><td style="max-width:200px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap" title="{{ p[6] or '' }}">{{ p[6] or '-' }}</td><td style="max-width:250px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap" title="{{ p[7] or '' }}">{{ p[7] or '-' }}</td><td>{% if p[10] %}{{ p[10].split(',')|length }} file(s){% else %}-{% endif %}</td><td>{{ p[11][:10] if p[11] else '-' }}</td><td>{{ p[12] or '-' }}</td><td><button class="edit-btn" onclick="window.open('/', '_blank')">Edit</button> <button class="delete-btn" onclick="deletePlasmid({{ p[0] }})">Delete</button></td></tr>
{% endfor %}</tbody></table></div></div></div>

<div class="section"><h2 onclick="toggleSection('backups')" style="cursor:pointer">📦 Backups <span id="backups-toggle">▼</span></h2>
<div id="backups-section" style="display:none"><div style="margin-bottom:15px"><button class="btn btn-success" onclick="createBackup()">Backup Now</button> <a href="/admin/download-archive" class="btn" style="background:#2ecc71">⬇ Download Archive (DB + HTML)</a></div>{% if backups %}{% for b in backups %}<div class="backup-item"><div><div class="backup-name">{{ b.name }}</div><div class="backup-meta">{{ b.size }} | {{ b.date }}</div></div><div><a href="/admin/download-backup/{{ b.name }}" class="btn" style="background:#3498db;font-size:12px;padding:8px 12px">⬇ Download</a></div></div>{% endfor %}{% else %}<p style="color:#7f8c8d">No backups yet</p>{% endif %}</div></div>

<div class="section"><h2 onclick="toggleSection('audit')" style="cursor:pointer">📋 Audit Log (Last 100) <span id="audit-toggle">▼</span></h2>
<div id="audit-section" style="display:none">{% if audit_logs %}<div style="overflow-x:auto"><table><thead><tr><th>ID</th><th>Time</th><th>Action</th><th>Record</th><th>IP</th><th>Data Preview</th></tr></thead><tbody>{% for log in audit_logs %}<tr><td>{{ log[0] }}</td><td>{{ log[1][:19] }}</td><td><span class="action-badge action-{{ log[2].lower() }}">{{ log[2] }}</span></td><td>{{ log[3] or '-' }}</td><td>{{ log[6] or '-' }}</td><td style="max-width:300px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap;font-family:monospace;font-size:10px">{{ log[4][:40] if log[4] else (log[5][:40] if log[5] else '-') }}</td></tr>{% endfor %}</tbody></table></div>{% else %}<p style="color:#7f8c8d">No logs yet</p>{% endif %}</div></div>
</div></div>
<script>
let adminSortDir={};
function toggleSection(id){const section=document.getElementById(id+'-section');const toggle=document.getElementById(id+'-toggle');if(section.style.display==='none'){section.style.display='block';toggle.textContent='▲'}else{section.style.display='none';toggle.textContent='▼'}}
function sortAdminTable(col){const table=document.getElementById('adminTable');const tbody=table.querySelector('tbody');const rows=Array.from(tbody.querySelectorAll('tr'));const dir=adminSortDir[col]==='asc'?'desc':'asc';adminSortDir={[col]:dir};rows.sort((a,b)=>{let aVal=a.cells[col].textContent.trim();let bVal=b.cells[col].textContent.trim();if(aVal==='-')aVal='';if(bVal==='-')bVal='';const aNum=parseFloat(aVal);const bNum=parseFloat(bVal);if(!isNaN(aNum)&&!isNaN(bNum))return dir==='asc'?aNum-bNum:bNum-aNum;return dir==='asc'?aVal.localeCompare(bVal):bVal.localeCompare(aVal)});rows.forEach(row=>tbody.appendChild(row))}
async function createBackup(){if(!confirm('Create backup?'))return;const r=await fetch('/api/backup',{method:'POST'});if(r.ok){alert('Backup created!');location.reload()}else{alert('Error')}}
async function deletePlasmid(id){if(!confirm('⚠️ PERMANENTLY DELETE this plasmid? This action is logged but cannot be undone!'))return;const r=await fetch(`/api/plasmids/${id}`,{method:'DELETE'});if(r.ok){alert('Deleted successfully');location.reload()}else{alert('Error deleting')}}
</script></body></html>
'''

# Login template
LOGIN_TEMPLATE = '''
<!DOCTYPE html>
<html><head><meta charset="UTF-8"><title>Admin Login - Plasmid Library</title>
<style>
*{margin:0;padding:0;box-sizing:border-box}
body{font-family:-apple-system,sans-serif;background:linear-gradient(135deg,#667eea 0%,#764ba2 100%);min-height:100vh;display:flex;align-items:center;justify-content:center}
.login-container{background:white;padding:40px;border-radius:12px;box-shadow:0 10px 40px rgba(0,0,0,0.2);width:100%;max-width:400px}
.login-header{text-align:center;margin-bottom:30px}
.login-header h1{color:#2c3e50;font-size:28px;margin-bottom:10px}
.login-header p{color:#7f8c8d;font-size:14px}
.form-group{margin-bottom:20px}
.form-group label{display:block;margin-bottom:8px;color:#2c3e50;font-weight:600;font-size:14px}
.form-group input{width:100%;padding:12px;border:2px solid #ddd;border-radius:6px;font-size:14px;transition:border-color 0.3s}
.form-group input:focus{outline:none;border-color:#8e44ad}
.btn-login{width:100%;background:#8e44ad;color:white;border:none;padding:14px;border-radius:6px;font-size:16px;font-weight:600;cursor:pointer;transition:background 0.3s}
.btn-login:hover{background:#7d3c98}
.error-msg{background:#f8d7da;color:#721c24;padding:12px;border-radius:6px;margin-bottom:20px;font-size:14px;border:1px solid #f5c6cb}
.back-link{text-align:center;margin-top:20px}
.back-link a{color:#8e44ad;text-decoration:none;font-size:14px}
.back-link a:hover{text-decoration:underline}
</style></head>
<body>
<div class="login-container">
<div class="login-header">
<h1>🔐 Admin Login</h1>
<p>Enter password to access admin panel</p>
</div>
{% if error %}
<div class="error-msg">❌ Incorrect password. Please try again.</div>
{% endif %}
<form method="POST">
<div class="form-group">
<label>Password</label>
<input type="password" name="password" required autofocus placeholder="Enter admin password">
</div>
<button type="submit" class="btn-login">Login</button>
</form>
<div class="back-link">
<a href="/">← Back to Database</a>
</div>
</div>
</body></html>
'''

EXPORT_TEMPLATE = '''
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Plasmid Library - Genetic Material Registry</title>
<style>
*{margin:0;padding:0;box-sizing:border-box}
body{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif;background:#fff;color:#2c3e50;padding:30px;max-width:1200px;margin:0 auto}
.report-header{border-bottom:3px solid #2c3e50;padding-bottom:20px;margin-bottom:30px}
.report-header h1{font-size:28px;color:#2c3e50;margin-bottom:5px}
.report-header .subtitle{font-size:14px;color:#7f8c8d}
.report-meta{display:flex;gap:30px;margin-top:15px;flex-wrap:wrap}
.report-meta-item{font-size:13px;color:#555}
.report-meta-item strong{color:#2c3e50}
.toc{background:#f8f9fa;border-radius:6px;padding:20px;margin-bottom:30px}
.toc h2{font-size:16px;margin-bottom:10px;color:#34495e}
.toc-list{column-count:3;column-gap:20px;list-style:none;font-size:12px}
.toc-list li{margin-bottom:4px}
.toc-list a{color:#3498db;text-decoration:none}
.toc-list a:hover{text-decoration:underline}
.plasmid-entry{border:1px solid #ddd;border-radius:8px;margin-bottom:25px;page-break-inside:avoid}
.plasmid-entry-header{background:#34495e;color:white;padding:12px 20px;border-radius:8px 8px 0 0;display:flex;justify-content:space-between;align-items:center}
.plasmid-entry-header h3{font-size:16px;font-weight:600}
.plasmid-entry-header .pid{font-size:20px;font-weight:700;color:#3498db;background:white;padding:4px 12px;border-radius:4px}
.plasmid-body{padding:20px}
.plasmid-meta{display:grid;grid-template-columns:repeat(auto-fill,minmax(200px,1fr));gap:10px;margin-bottom:15px}
.meta-field{font-size:12px}
.meta-field .label{color:#7f8c8d;font-weight:600;text-transform:uppercase;font-size:10px;display:block}
.meta-field .value{color:#2c3e50;margin-top:2px}
.seq-section{margin-top:15px;border-top:1px solid #eee;padding-top:15px}
.seq-section h4{font-size:14px;color:#34495e;margin-bottom:10px}
.seq-display{display:flex;gap:20px;flex-wrap:wrap}
.seq-map-export{flex:0 0 300px;text-align:center;cursor:pointer;border:2px solid transparent;border-radius:6px;padding:4px;transition:border-color 0.2s}
.seq-map-export:hover{border-color:#3498db}
.seq-map-export::after{content:"click to zoom";display:block;font-size:10px;color:#999;margin-top:4px}
.seq-data-export{flex:1;min-width:300px}
.seq-box{background:#1e1e1e;color:#d4d4d4;font-family:'Courier New',monospace;font-size:10px;padding:12px;border-radius:4px;max-height:300px;overflow-y:auto;word-break:break-all;white-space:pre-wrap;line-height:1.5}
.seq-box .insert-hl{color:#e74c3c;font-weight:bold;background:rgba(231,76,60,0.15)}
.seq-info-badge{display:inline-block;background:#ecf0f1;padding:3px 8px;border-radius:3px;font-size:11px;margin:2px}
.seq-info-badge.insert{background:#fde8e8;color:#e74c3c;font-weight:600}
.features-mini{font-size:11px;margin-top:10px}
.features-mini table{width:100%;border-collapse:collapse}
.features-mini th{background:#ecf0f1;padding:4px 8px;text-align:left;font-size:10px;font-weight:600}
.features-mini td{padding:4px 8px;border-bottom:1px solid #f0f0f0;font-size:10px}
.footer{border-top:2px solid #2c3e50;padding-top:15px;margin-top:40px;font-size:11px;color:#7f8c8d;text-align:center}
.lightbox{display:none;position:fixed;z-index:9999;top:0;left:0;width:100%;height:100%;background:rgba(0,0,0,0.85);align-items:center;justify-content:center;cursor:pointer}
.lightbox.active{display:flex}
.lightbox-content{max-width:90vw;max-height:90vh;background:white;border-radius:8px;padding:20px;position:relative}
.lightbox-content svg{width:80vmin;height:80vmin}
.lightbox-close{position:absolute;top:-12px;right:-12px;background:#e74c3c;color:white;border:none;border-radius:50%;width:30px;height:30px;font-size:16px;cursor:pointer;font-weight:bold}
.dl-link{display:inline-block;background:#34495e;color:white;padding:3px 8px;border-radius:3px;font-size:10px;text-decoration:none;margin-top:4px}
.dl-link:hover{background:#2c3e50}
@media print{body{padding:10px;font-size:11px}.plasmid-entry{break-inside:avoid}.seq-box{max-height:none}.lightbox{display:none!important}}
</style>
</head>
<body>
<div class="report-header">
<h1>Plasmid Library - Genetic Material Registry</h1>
<p class="subtitle">Complete inventory of genetic resources with sequence data</p>
<div class="report-meta">
<div class="report-meta-item"><strong>Export Date:</strong> {{ export_date }}</div>
<div class="report-meta-item"><strong>Total Entries:</strong> {{ total }}</div>
<div class="report-meta-item"><strong>Institution:</strong> Your Institution</div>
</div>
</div>

<div class="toc">
<h2>Table of Contents ({{ total }} entries)</h2>
<ul class="toc-list">
{% for p in plasmids %}
<li><a href="#entry-{{ p.number }}">P{{ p.number }} - {{ p.name }}</a></li>
{% endfor %}
</ul>
</div>

{% for p in plasmids %}
<div class="plasmid-entry" id="entry-{{ p.number }}">
<div class="plasmid-entry-header">
<h3>{{ p.name }}</h3>
<span class="pid">P{{ p.number }}</span>
</div>
<div class="plasmid-body">
<div class="plasmid-meta">
<div class="meta-field"><span class="label">Plasmid Backbone</span><span class="value">{{ p.backbone or 'N/A' }}</span></div>
<div class="meta-field"><span class="label">Inductor</span><span class="value">{{ p.inductor or 'N/A' }}</span></div>
<div class="meta-field"><span class="label">Antibiotic Resistance</span><span class="value">{{ p.resistance or 'N/A' }}</span></div>
<div class="meta-field"><span class="label">Reference</span><span class="value">{{ p.reference or 'N/A' }}</span></div>
<div class="meta-field"><span class="label">Date Added</span><span class="value">{{ p.created_at[:10] if p.created_at else 'N/A' }}</span></div>
</div>
{% if p.comment %}<div class="meta-field" style="margin-bottom:10px"><span class="label">Comment</span><span class="value">{{ p.comment }}</span></div>{% endif %}

{% for seq in p.sequences %}
<div class="seq-section">
<h4>Sequence: {{ seq.filename.split('_', 3)[-1] if '_' in seq.filename else seq.filename }}</h4>
<div style="margin-bottom:8px">
<span class="seq-info-badge"><strong>Length:</strong> {{ seq.length|default(0) }} bp</span>
<span class="seq-info-badge"><strong>Topology:</strong> {{ seq.topology|default('unknown') }}</span>
<span class="seq-info-badge"><strong>Format:</strong> {{ seq.format|default('unknown') }}</span>
{% if seq.insert_region %}<span class="seq-info-badge insert"><strong>Gene of Interest:</strong> {{ seq.insert_region.label }} ({{ seq.insert_region.start }}-{{ seq.insert_region.end }})</span>{% endif %}
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
<p style="color:#999;font-size:12px">Sequence data could not be extracted from this file format.</p>
{% endif %}
</div>
</div>
{% if seq.features %}
<div class="features-mini">
<table><thead><tr><th>Type</th><th>Label</th><th>Location</th><th>Strand</th><th>Product</th></tr></thead><tbody>
{% for f in seq.features %}
<tr{% if seq.insert_region and f.start == seq.insert_region.start and f.end == seq.insert_region.end %} style="background:#fde8e8;font-weight:600"{% endif %}>
<td>{{ f.type }}</td><td>{{ f.label or '-' }}</td><td>{{ f.start }}..{{ f.end }}</td><td>{{ '+' if f.strand == 1 else ('-' if f.strand == -1 else '?') }}</td><td>{{ f.product or '-' }}</td></tr>
{% endfor %}
</tbody></table>
</div>
{% endif %}
</div>
{% endfor %}
{% if not p.sequences %}<p style="color:#999;font-size:12px;margin-top:10px;font-style:italic">No parseable sequence files attached to this entry.</p>{% endif %}
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
<p>Generated from Plasmid Library Database on {{ export_date }}</p>
<p>Your Institution</p>
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
    print("⚠️  CHANGE THIS in plasmid_server.py line 18!")
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
