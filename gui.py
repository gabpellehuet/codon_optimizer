"""
Codon Optimizer - GUI
======================
Tkinter front-end for reverse_translation.py. Pure stdlib (tkinter ships with
Python), so it can be frozen into a single Windows executable with PyInstaller
without pulling in any extra dependencies.

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

CUSTOM_TABLE_LABEL = "Custom / loaded from file"

DEFAULT_AA_SEQUENCE = (
    "EIIEEVLLNSKKDSLTEAEIYELVQKQLSTDPALKGVKLSKEEVRETLQELVIEGRLIKDKITGKYRLSTNTRLELLIEQLP"
)


class CodonOptimizerApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Codon Optimizer")
        self.geometry("1000x720")
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
        organism_names = list(CODON_TABLES.keys()) + [CUSTOM_TABLE_LABEL]
        self.organism_combo = ttk.Combobox(
            picker_row, textvariable=self.organism, values=organism_names,
            state="readonly", width=28,
        )
        self.organism_combo.pack(side="left", padx=(5, 0))
        self.organism_combo.bind("<<ComboboxSelected>>", self._on_organism_change)

        ttk.Button(table_frame, text="Load from file...",
                   command=self._load_table_file).grid(row=1, column=0, sticky="w", pady=(5, 0))
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

    def _on_mode_change(self):
        self.seq_text.delete("1.0", "end")
        if self.mode.get() == "aa":
            self.seq_text.insert("1.0", DEFAULT_AA_SEQUENCE)

    def _on_organism_change(self, event=None):
        name = self.organism.get()
        if name == CUSTOM_TABLE_LABEL:
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
