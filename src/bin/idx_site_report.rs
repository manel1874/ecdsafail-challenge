use std::collections::BTreeMap;
use std::io::Write;

use quantum_ecc::circuit::OperationType;
use quantum_ecc::point_add;

fn main() {
    std::env::set_var("TRACE_OP_SITES", "1");
    let idx_path = std::env::var("IDX_SITE_FILE").expect("set IDX_SITE_FILE=/path/to.idx");
    let include_context = std::env::var("IDX_SITE_CONTEXT").ok().as_deref() == Some("1");
    let split_dir = std::env::var("IDX_SITE_SPLIT_DIR").ok();
    let limit = std::env::var("IDX_SITE_LIMIT")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(80);
    let indices = read_indices(&idx_path);
    let ops = point_add::build();
    let sites = point_add::take_last_op_sites();
    let mut counts = BTreeMap::<String, usize>::new();
    let mut samples = BTreeMap::<String, Vec<usize>>::new();
    for idx in indices {
        if idx >= ops.len() {
            eprintln!("idx_site_report: skipping out-of-range idx={idx}");
            continue;
        }
        if !matches!(ops[idx].kind, OperationType::CCX | OperationType::CCZ) {
            eprintln!(
                "idx_site_report: idx={} is {:?}, not CCX/CCZ",
                idx, ops[idx].kind
            );
        }
        let site = format_site(&sites, idx, include_context);
        *counts.entry(site.clone()).or_insert(0) += 1;
        let sample = samples.entry(site).or_default();
        if sample.len() < 12 {
            sample.push(idx);
        }
    }
    if let Some(dir) = split_dir.as_ref() {
        std::fs::create_dir_all(dir)
            .unwrap_or_else(|err| panic!("create IDX_SITE_SPLIT_DIR={dir}: {err}"));
        for (site, items) in &samples {
            let _ = (site, items);
        }
        let mut grouped = BTreeMap::<String, Vec<usize>>::new();
        let indices = read_indices(&idx_path);
        for idx in indices {
            if idx < ops.len() {
                grouped
                    .entry(format_site(&sites, idx, include_context))
                    .or_default()
                    .push(idx);
            }
        }
        for (site, items) in grouped {
            let path = format!("{}/{}.idx", dir, sanitize_filename(&site));
            let mut f = std::fs::File::create(&path)
                .unwrap_or_else(|err| panic!("create split idx {path}: {err}"));
            writeln!(f, "idx").unwrap();
            for idx in items {
                writeln!(f, "{idx}").unwrap();
            }
            f.flush().unwrap();
        }
    }
    let mut rows: Vec<_> = counts.into_iter().collect();
    rows.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
    println!(
        "idx_site_report idx_file={} groups={} showing_top={}",
        idx_path,
        rows.len(),
        limit.min(rows.len())
    );
    for (site, count) in rows.into_iter().take(limit) {
        let sample = samples
            .get(&site)
            .map(|items| {
                items
                    .iter()
                    .map(|idx| idx.to_string())
                    .collect::<Vec<_>>()
                    .join(",")
            })
            .unwrap_or_default();
        println!("idx_site count={count} site={site} sample={sample}");
    }
}

fn read_indices(path: &str) -> Vec<usize> {
    let text = std::fs::read_to_string(path)
        .unwrap_or_else(|err| panic!("failed to read IDX_SITE_FILE={path}: {err}"));
    text.lines()
        .filter_map(|line| {
            let t = line.trim();
            if t.is_empty() || t.starts_with('#') || t.starts_with("idx") {
                return None;
            }
            t.split('\t').next().unwrap_or(t).parse::<usize>().ok()
        })
        .collect()
}

fn format_site(sites: &[point_add::OpSite], idx: usize, include_context: bool) -> String {
    if sites.len() <= idx {
        return "unavailable".to_owned();
    }
    let (file, line, context) = sites[idx];
    if include_context && context != 0 {
        format!("{file}:{line}#0x{context:08x}")
    } else {
        format!("{file}:{line}")
    }
}

fn sanitize_filename(value: &str) -> String {
    value
        .chars()
        .map(|c| match c {
            'A'..='Z' | 'a'..='z' | '0'..='9' | '-' | '_' | '.' => c,
            _ => '_',
        })
        .collect()
}
