# Phase C review

Five lanes reviewed source/evidence. Goal reviewer found overlapping split/merge
suppression; code reviewer found primary inversion tp:A:I incorrectly skipped.
Both were reproduced (2978175: 18 pass, 4 fail), fixed, and independently
re-reviewed PASS. Security/context PASS; targeted QA PASS22/22 in2978176.
Four further boundary tests bring latest2978177 to26/26 without production edits.
Original legacy same-name single-indel comparator emitted4540 events vs1;
reverse-complement comparison emitted2 unmatched contig events. Actual minimap2
smoke tests now produce1 deletion and1 RC event respectively.
Dotted command labels also have reproduced/fixed receipt-collision coverage.
No Phase D whole-benchmark result is implied. Early smoke tests use temporary
receipts; latest source hashes/minimap2 version are logged, and formal execution
preserves command/PAF receipts permanently.
