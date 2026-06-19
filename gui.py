"""
Codon Optimizer - GUI
======================
Tkinter front-end for reverse_translation.py. The GUI itself only needs
tkinter (stdlib), but the restriction-site tab depends on Biopython
(Bio.Restriction) -- run `pip install biopython` before running from source
or rebuilding the frozen executable.

Run with:
    python gui.py
"""

from __future__ import annotations
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

from reverse_translation import (
    parse_codon_table,
    optimize_from_aa,
    optimize_from_nt,
    _aa_for_codon,
    CODON_TABLE_TEXT,
)
from codon_tables import CODON_TABLES
import kazusa_online
import restriction_sites

CUSTOM_TABLE_LABEL = "Custom / loaded from file"
ONLINE_TABLE_LABEL = "Custom / fetched from Kazusa online"

DEFAULT_AA_SEQUENCE = (
    "EIIEEVLLNSKKDSLTEAEIYELVQKQLSTDPALKGVKLSKEEVRETLQELVIEGRLIKDKITGKYRLSTNTRLELLIEQLP"
)


def _get_work_area():
    """Return (top, usable_height) of the screen area excluding the Windows
    taskbar, via the SPI_GETWORKAREA system call. Falls back to a rough
    estimate on non-Windows platforms (e.g. running the GUI from source on
    Linux/macOS for testing)."""
    try:
        import ctypes

        class RECT(ctypes.Structure):
            _fields_ = [("left", ctypes.c_long), ("top", ctypes.c_long),
                        ("right", ctypes.c_long), ("bottom", ctypes.c_long)]

        rect = RECT()
        SPI_GETWORKAREA = 0x0030
        ctypes.windll.user32.SystemParametersInfoW(SPI_GETWORKAREA, 0, ctypes.byref(rect), 0)
        return rect.top, rect.bottom - rect.top
    except (AttributeError, OSError):
        return 0, None


class CodonOptimizerApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Codon Optimizer")
        top, usable_h = _get_work_area()
        if usable_h is None:
            usable_h = self.winfo_screenheight() - 80
        TITLE_BAR_H = 32  # tkinter's geometry height excludes the title bar
        height = max(600, min(680, usable_h - TITLE_BAR_H - 20))
        self.geometry(f"1000x{height}+100+{top + 10}")
        self.minsize(800, 600)
        self._result = None
        self._build_widgets()

    # ── UI construction ────────────────────────────────────────────────
    def _build_widgets(self):
        # Top-level layout uses pack(), not grid(), so the Optimize button and
        # the bottom action bar are pinned at their natural size and never get
        # squeezed out when the window is resized smaller. Only the results
        # notebook (packed last, with expand=True) absorbs extra/missing space.
        root = ttk.Frame(self, padding=10)
        root.pack(fill="both", expand=True)

        # Mode selector (pinned, top)
        mode_frame = ttk.LabelFrame(root, text="Input mode", padding=8)
        mode_frame.pack(side="top", fill="x", pady=(0, 8))
        self.mode = tk.StringVar(value="aa")
        ttk.Radiobutton(mode_frame, text="Amino acid sequence", variable=self.mode,
                         value="aa", command=self._on_mode_change).pack(side="left", padx=5)
        ttk.Radiobutton(mode_frame, text="Nucleotide sequence", variable=self.mode,
                         value="nt", command=self._on_mode_change).pack(side="left", padx=5)

        ttk.Label(mode_frame, text="Frequency threshold:").pack(side="left", padx=(20, 5))
        self.freq_threshold = tk.DoubleVar(value=0.10)
        ttk.Spinbox(mode_frame, from_=0.0, to=1.0, increment=0.01, width=6,
                    textvariable=self.freq_threshold).pack(side="left")

        # Input row: sequence + codon table side by side (pinned, top)
        input_frame = ttk.Frame(root)
        input_frame.pack(side="top", fill="x", pady=(0, 8))
        input_frame.columnconfigure(0, weight=1)
        input_frame.columnconfigure(1, weight=1)

        seq_frame = ttk.LabelFrame(input_frame, text="Sequence", padding=8)
        seq_frame.grid(row=0, column=0, sticky="nsew", padx=(0, 5))
        self.seq_text = tk.Text(seq_frame, height=5, wrap="word")
        self.seq_text.pack(fill="both", expand=True)
        self.seq_text.insert("1.0", DEFAULT_AA_SEQUENCE)

        table_frame = ttk.LabelFrame(input_frame, text="Codon usage table", padding=8)
        table_frame.grid(row=0, column=1, sticky="nsew", padx=(5, 0))
        table_frame.columnconfigure(0, weight=1)
        table_frame.rowconfigure(2, weight=1)

        picker_row = ttk.Frame(table_frame)
        picker_row.grid(row=0, column=0, sticky="ew")
        ttk.Label(picker_row, text="Organism:").pack(side="left")
        self.organism = tk.StringVar(value="Escherichia coli K12")
        organism_names = list(CODON_TABLES.keys()) + [CUSTOM_TABLE_LABEL, ONLINE_TABLE_LABEL]
        self.organism_combo = ttk.Combobox(
            picker_row, textvariable=self.organism, values=organism_names,
            state="readonly", width=28,
        )
        self.organism_combo.pack(side="left", padx=(5, 0))
        self.organism_combo.bind("<<ComboboxSelected>>", self._on_organism_change)

        button_row = ttk.Frame(table_frame)
        button_row.grid(row=1, column=0, sticky="w", pady=(5, 0))
        ttk.Button(button_row, text="Load from file...",
                   command=self._load_table_file).pack(side="left")
        ttk.Button(button_row, text="Fetch from Kazusa online...",
                   command=self._fetch_from_kazusa).pack(side="left", padx=(5, 0))

        self.table_text = tk.Text(table_frame, height=5, wrap="none")
        self.table_text.grid(row=2, column=0, sticky="nsew", pady=(5, 0))
        self.table_text.insert("1.0", CODON_TABLE_TEXT.strip())

        # Bottom action bar (pinned, bottom) - packed before the expanding
        # notebook so it always keeps its place at the very bottom.
        action_frame = ttk.Frame(root)
        action_frame.pack(side="bottom", fill="x", pady=(8, 0))
        ttk.Button(action_frame, text="Copy final sequence",
                   command=self._copy_sequence).pack(side="left")
        ttk.Button(action_frame, text="Save to FASTA...",
                   command=self._save_fasta).pack(side="left", padx=5)

        # Run button (pinned, just above the action bar)
        run_btn = ttk.Button(root, text="Optimize", command=self._run_optimization)
        run_btn.pack(side="bottom", fill="x", pady=(0, 8))

        # Results notebook - the only section that expands/shrinks
        self.notebook = ttk.Notebook(root)
        self.notebook.pack(side="top", fill="both", expand=True)

        self.summary_text = tk.Text(self.notebook, height=8, wrap="word", state="disabled")
        self.notebook.add(self.summary_text, text="Sequences")

        self.changes_tree = ttk.Treeview(
            self.notebook, columns=("idx", "aa", "old", "new", "reason"), show="headings"
        )
        for col, label, width in [
            ("idx", "Codon #", 70), ("aa", "AA", 60), ("old", "Original", 80),
            ("new", "Replaced", 80), ("reason", "Reason", 400),
        ]:
            self.changes_tree.heading(col, text=label)
            self.changes_tree.column(col, width=width, anchor="center")
        self.notebook.add(self.changes_tree, text="Changes")

        self.issues_text = tk.Text(self.notebook, height=8, wrap="word", state="disabled")
        self.notebook.add(self.issues_text, text="Warnings / Remaining issues")

        restriction_tab = ttk.Frame(self.notebook)
        self.notebook.add(restriction_tab, text="Restriction Sites")
        self._build_restriction_tab(restriction_tab)

    def _build_restriction_tab(self, parent):
        # Available enzymes (filterable), left
        avail_frame = ttk.Frame(parent, padding=8)
        avail_frame.pack(side="left", fill="y")
        ttk.Label(avail_frame, text="Available enzymes").pack(side="top", anchor="w")
        ttk.Label(avail_frame, text="Filter:").pack(side="top", anchor="w")
        self.enzyme_filter = tk.StringVar()
        filter_entry = ttk.Entry(avail_frame, textvariable=self.enzyme_filter, width=20)
        filter_entry.pack(side="top", fill="x")
        filter_entry.bind("<KeyRelease>", self._on_enzyme_filter_change)

        avail_list_frame = ttk.Frame(avail_frame)
        avail_list_frame.pack(side="top", fill="both", expand=True, pady=(5, 0))
        self.enzyme_listbox = tk.Listbox(avail_list_frame, selectmode="extended",
                                          height=8, exportselection=False)
        self.enzyme_listbox.pack(side="left", fill="both", expand=True)
        avail_scroll = ttk.Scrollbar(avail_list_frame, command=self.enzyme_listbox.yview)
        avail_scroll.pack(side="right", fill="y")
        self.enzyme_listbox.config(yscrollcommand=avail_scroll.set)
        self.enzyme_listbox.bind("<Double-Button-1>", lambda e: self._add_selected_enzymes())

        self._all_enzyme_names = restriction_sites.list_enzyme_names()
        self._populate_enzyme_listbox(self._all_enzyme_names)

        # Add / Remove buttons, middle
        move_frame = ttk.Frame(parent, padding=8)
        move_frame.pack(side="left", fill="y")
        ttk.Button(move_frame, text="Add →",
                   command=self._add_selected_enzymes).pack(side="top", pady=(60, 5))
        ttk.Button(move_frame, text="← Remove",
                   command=self._remove_selected_enzymes).pack(side="top")

        # Selected enzymes + action buttons, next
        sel_frame = ttk.Frame(parent, padding=8)
        sel_frame.pack(side="left", fill="y")
        ttk.Label(sel_frame, text="Selected enzymes").pack(side="top", anchor="w")
        sel_list_frame = ttk.Frame(sel_frame)
        sel_list_frame.pack(side="top", fill="both", expand=True, pady=(5, 0))
        self.selected_listbox = tk.Listbox(sel_list_frame, selectmode="extended",
                                            height=8, exportselection=False)
        self.selected_listbox.pack(side="left", fill="both", expand=True)
        sel_scroll = ttk.Scrollbar(sel_list_frame, command=self.selected_listbox.yview)
        sel_scroll.pack(side="right", fill="y")
        self.selected_listbox.config(yscrollcommand=sel_scroll.set)
        self.selected_listbox.bind("<Double-Button-1>", lambda e: self._remove_selected_enzymes())

        button_col = ttk.Frame(sel_frame)
        button_col.pack(side="top", fill="x", pady=(5, 0))
        ttk.Button(button_col, text="Show restriction sites",
                   command=self._show_restriction_sites).pack(side="top", fill="x")
        ttk.Button(button_col, text="Eliminate restriction sites",
                   command=self._eliminate_restriction_sites).pack(side="top", fill="x", pady=(3, 0))

        # Results, right (expands)
        right = ttk.Frame(parent, padding=8)
        right.pack(side="left", fill="both", expand=True)
        self.restriction_tree = ttk.Treeview(
            right, columns=("enzyme", "position", "site"), show="headings"
        )
        for col, label, width in [("enzyme", "Enzyme", 100),
                                   ("position", "Position (nt)", 100),
                                   ("site", "Recognition site", 150)]:
            self.restriction_tree.heading(col, text=label)
            self.restriction_tree.column(col, width=width, anchor="center")
        self.restriction_tree.pack(fill="both", expand=True)

    def _populate_enzyme_listbox(self, names):
        self.enzyme_listbox.delete(0, "end")
        for name in names:
            self.enzyme_listbox.insert("end", name)

    def _on_enzyme_filter_change(self, event=None):
        query = self.enzyme_filter.get().strip().lower()
        names = [n for n in self._all_enzyme_names if query in n.lower()] if query \
            else self._all_enzyme_names
        self._populate_enzyme_listbox(names)

    def _add_selected_enzymes(self):
        chosen = [self.enzyme_listbox.get(i) for i in self.enzyme_listbox.curselection()]
        current = set(self.selected_listbox.get(0, "end"))
        current.update(chosen)
        self.selected_listbox.delete(0, "end")
        for name in sorted(current):
            self.selected_listbox.insert("end", name)

    def _remove_selected_enzymes(self):
        sel = list(self.selected_listbox.curselection())
        for i in reversed(sel):
            self.selected_listbox.delete(i)

    def _selected_enzymes(self):
        return list(self.selected_listbox.get(0, "end"))

    def _refresh_restriction_tree(self, sites: dict):
        for row in self.restriction_tree.get_children():
            self.restriction_tree.delete(row)
        final_nt = "".join(self._result["final_codons"]) if self._result else ""
        for enzyme_name, spans in sites.items():
            for start, end in spans:
                self.restriction_tree.insert(
                    "", "end", values=(enzyme_name, start, final_nt[start:end])
                )

    def _show_restriction_sites(self):
        if not self._result:
            messagebox.showinfo("Nothing to scan", "Run an optimization first.")
            return
        enzymes = self._selected_enzymes()
        if not enzymes:
            messagebox.showinfo("No enzymes selected",
                                 "Select one or more restriction enzymes from the list first.")
            return
        sites = restriction_sites.find_sites(self._result["final_codons"], enzymes)
        self._refresh_restriction_tree(sites)
        self.notebook.select(3)

    def _eliminate_restriction_sites(self):
        if not self._result:
            messagebox.showinfo("Nothing to fix", "Run an optimization first.")
            return
        enzymes = self._selected_enzymes()
        if not enzymes:
            messagebox.showinfo("No enzymes selected",
                                 "Select one or more restriction enzymes from the list first.")
            return
        ct = self._result["codon_table"]
        threshold = self.freq_threshold.get()
        new_codons, changes, warnings, remaining = restriction_sites.eliminate_sites(
            self._result["final_codons"], ct, enzymes, freq_threshold=threshold,
        )
        self._result["final_codons"] = new_codons
        self._result["changes"] = self._result["changes"] + changes
        self._result["warnings"] = self._result["warnings"] + warnings
        self._display_result(self._result, threshold)
        self._refresh_restriction_tree(remaining)
        self.notebook.select(3)

    def _on_mode_change(self):
        self.seq_text.delete("1.0", "end")
        if self.mode.get() == "aa":
            self.seq_text.insert("1.0", DEFAULT_AA_SEQUENCE)

    def _on_organism_change(self, event=None):
        name = self.organism.get()
        if name in (CUSTOM_TABLE_LABEL, ONLINE_TABLE_LABEL):
            return
        self.table_text.delete("1.0", "end")
        self.table_text.insert("1.0", CODON_TABLES[name].strip())

    def _load_table_file(self):
        path = filedialog.askopenfilename(
            title="Select codon usage table file",
            filetypes=[("Text files", "*.txt"), ("All files", "*.*")],
        )
        if not path:
            return
        with open(path, "r", encoding="utf-8") as f:
            content = f.read()
        self.table_text.delete("1.0", "end")
        self.table_text.insert("1.0", content)
        self.organism.set(CUSTOM_TABLE_LABEL)

    def _fetch_from_kazusa(self):
        dialog = tk.Toplevel(self)
        dialog.title("Fetch from Kazusa Codon Usage Database")
        dialog.geometry("520x400")
        dialog.transient(self)

        top_row = ttk.Frame(dialog, padding=8)
        top_row.pack(side="top", fill="x")
        ttk.Label(top_row, text="Organism name:").pack(side="left")
        query_var = tk.StringVar()
        entry = ttk.Entry(top_row, textvariable=query_var, width=30)
        entry.pack(side="left", padx=5, fill="x", expand=True)
        entry.focus_set()

        status_var = tk.StringVar(value="Type an organism name and press Search.")
        status_label = ttk.Label(dialog, textvariable=status_var, padding=(8, 0))
        status_label.pack(side="top", fill="x")

        results_frame = ttk.Frame(dialog, padding=8)
        results_frame.pack(side="top", fill="both", expand=True)
        listbox = tk.Listbox(results_frame)
        listbox.pack(side="left", fill="both", expand=True)
        scrollbar = ttk.Scrollbar(results_frame, command=listbox.yview)
        scrollbar.pack(side="right", fill="y")
        listbox.config(yscrollcommand=scrollbar.set)

        button_row = ttk.Frame(dialog, padding=8)
        button_row.pack(side="bottom", fill="x")

        search_results: list[tuple[str, str]] = []

        def do_search(event=None):
            query = query_var.get().strip()
            if not query:
                return
            status_var.set("Searching...")
            dialog.update()
            try:
                results = kazusa_online.search_organisms(query)
            except Exception as exc:
                status_var.set(f"Search failed: {exc}")
                return
            search_results.clear()
            search_results.extend(results)
            listbox.delete(0, "end")
            for _, label in results:
                listbox.insert("end", label)
            status_var.set(f"{len(results)} result(s) found." if results
                            else "No matches found.")

        def do_use_selected():
            sel = listbox.curselection()
            if not sel:
                messagebox.showinfo("No selection", "Select an organism from the list first.")
                return
            species_id, label = search_results[sel[0]]
            status_var.set(f"Fetching codon table for {label}...")
            dialog.update()
            try:
                table_text = kazusa_online.fetch_codon_table(species_id)
            except Exception as exc:
                messagebox.showerror("Fetch failed", str(exc))
                return
            self.table_text.delete("1.0", "end")
            self.table_text.insert("1.0", table_text)
            self.organism.set(ONLINE_TABLE_LABEL)
            dialog.destroy()

        entry.bind("<Return>", do_search)
        ttk.Button(top_row, text="Search", command=do_search).pack(side="left")
        ttk.Button(button_row, text="Use selected",
                   command=do_use_selected).pack(side="left")
        ttk.Button(button_row, text="Cancel",
                   command=dialog.destroy).pack(side="left", padx=5)

    # ── Core action ─────────────────────────────────────────────────────
    def _run_optimization(self):
        seq = self.seq_text.get("1.0", "end").strip()
        table_text = self.table_text.get("1.0", "end")
        threshold = self.freq_threshold.get()

        if not seq:
            messagebox.showerror("Missing input", "Please enter a sequence.")
            return

        try:
            parse_codon_table(table_text)  # validate early
            if self.mode.get() == "aa":
                result = optimize_from_aa(seq, table_text, freq_threshold=threshold)
            else:
                result = optimize_from_nt(seq, table_text, freq_threshold=threshold)
        except Exception as exc:
            messagebox.showerror("Optimization failed", str(exc))
            return

        self._result = result
        self._display_result(result, threshold)

    # ── Result display ─────────────────────────────────────────────────
    def _display_result(self, result: dict, threshold: float):
        ct = result["codon_table"]
        initial = result["initial_codons"]
        final = result["final_codons"]
        changes = result["changes"]
        warnings = result["warnings"]
        issues = result["remaining_issues"]
        initial_nt = "".join(initial)
        final_nt = "".join(final)

        # Summary tab
        self.summary_text.config(state="normal")
        self.summary_text.delete("1.0", "end")
        self.summary_text.insert("end", f"Length: {len(final_nt)} nt ({len(final)} codons)\n")
        self.summary_text.insert("end", f"Changes applied: {len(changes)}\n")
        self.summary_text.insert("end", f"Warnings: {len(warnings)}\n")
        self.summary_text.insert("end", f"Frequency threshold: {threshold}\n\n")
        self.summary_text.insert("end", "Initial sequence (before optimization):\n")
        self.summary_text.insert("end", initial_nt + "\n\n")
        self.summary_text.insert("end", "Final sequence (after optimization):\n")
        self.summary_text.insert("end", final_nt + "\n")
        self.summary_text.config(state="disabled")

        # Changes tab
        for row in self.changes_tree.get_children():
            self.changes_tree.delete(row)
        for idx, old, new, reason in changes:
            aa = _aa_for_codon(old, ct)
            self.changes_tree.insert("", "end", values=(idx, aa, old, new, reason))

        # Issues/warnings tab
        self.issues_text.config(state="normal")
        self.issues_text.delete("1.0", "end")
        has_issues = any(len(v) > 0 for v in issues.values())
        if has_issues:
            self.issues_text.insert("end", "REMAINING ISSUES:\n")
            for i, j, c in issues.get("identical_adjacent", []):
                aa = _aa_for_codon(c, ct)
                self.issues_text.insert("end", f"  Adjacent: codons {i} & {j}: {c}|{c} ({aa})\n")
            for i, j, c in issues.get("identical_skip1", []):
                mid = final[i + 1]
                aa = _aa_for_codon(c, ct)
                self.issues_text.insert("end", f"  Skip-1: codons {i} & {j}: {c}_{mid}_{c} ({aa})\n")
            for nt, start, run_len in issues.get("homopolymer_runs", []):
                self.issues_text.insert("end", f"  Homopolymer: '{nt * run_len}' @{start}\n")
        else:
            self.issues_text.insert("end", "No remaining issues.\n")
        if warnings:
            self.issues_text.insert("end", "\nWARNINGS:\n")
            for w in warnings:
                self.issues_text.insert("end", f"  - {w}\n")
        self.issues_text.config(state="disabled")

    # ── Output actions ─────────────────────────────────────────────────
    def _copy_sequence(self):
        if not self._result:
            messagebox.showinfo("Nothing to copy", "Run an optimization first.")
            return
        final_nt = "".join(self._result["final_codons"])
        self.clipboard_clear()
        self.clipboard_append(final_nt)
        messagebox.showinfo("Copied", "Final sequence copied to clipboard.")

    def _save_fasta(self):
        if not self._result:
            messagebox.showinfo("Nothing to save", "Run an optimization first.")
            return
        path = filedialog.asksaveasfilename(
            title="Save FASTA", defaultextension=".fasta",
            filetypes=[("FASTA files", "*.fasta"), ("All files", "*.*")],
        )
        if not path:
            return
        final_nt = "".join(self._result["final_codons"])
        with open(path, "w", encoding="utf-8") as f:
            f.write(">optimized_sequence\n")
            for i in range(0, len(final_nt), 70):
                f.write(final_nt[i:i + 70] + "\n")
        messagebox.showinfo("Saved", f"Saved to {path}")


def main():
    app = CodonOptimizerApp()
    app.mainloop()


if __name__ == "__main__":
    main()
