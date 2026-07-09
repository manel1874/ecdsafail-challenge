use quantum_ecc::circuit::{Op, NO_BIT, NO_QUBIT};
use quantum_ecc::point_add;

fn main() {
    let start = env_usize("OP_WINDOW_START").unwrap_or(0);
    let len = env_usize("OP_WINDOW_LEN").unwrap_or(32);
    let ops = point_add::build();
    let end = start.saturating_add(len).min(ops.len());
    println!("op_window ops={} range={}..{}", ops.len(), start, end);
    for (idx, op) in ops.iter().enumerate().take(end).skip(start) {
        println!("{idx}: {}", format_op(*op));
    }
}

fn env_usize(name: &str) -> Option<usize> {
    std::env::var(name)
        .ok()
        .and_then(|value| value.parse::<usize>().ok())
}

fn format_op(op: Op) -> String {
    let mut parts = vec![format!("{:?}", op.kind)];
    for (label, q) in [
        ("c2", op.q_control2),
        ("c1", op.q_control1),
        ("t", op.q_target),
    ] {
        if q != NO_QUBIT {
            parts.push(format!("{label}=q{}", q.0));
        }
    }
    if op.c_target != NO_BIT {
        parts.push(format!("bt=b{}", op.c_target.0));
    }
    if op.c_condition != NO_BIT {
        parts.push(format!("if=b{}", op.c_condition.0));
    }
    parts.join(" ")
}
