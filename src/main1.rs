fn mod_inv(a: i64, p: i64) -> i64 {
    let (mut old_r, mut r) = (a, p);
    let (mut old_t, mut t) = (1, 0);

    while r != 0 {
        let quotient = old_r / r;
        old_r = old_r - quotient * r;
        std::mem::swap(&mut old_r, &mut r);

        old_t = old_t - quotient * t;
        std::mem::swap(&mut old_t, &mut t);
    }

    if old_t < 0 {
        old_t += p;
    }

    old_t
}

use rand::Rng;

fn generate_shares(
    secret: i64,
    threshold: usize,
    num_shares: usize,
    prime: i64,
) -> Vec<(i64, i64)> {
    let mut rng = rand::rng();
    println!("range: {:?}", rng.clone());
    // Random coefficients (except a_0 which is the secret)
    let mut coefficients = vec![secret];
    for x in 1..threshold {
        coefficients.push(rng.random_range(1..prime));
        println!("coefficient {} is {:?}", x, coefficients[x]);
    }

    // Generate shares
    let mut shares = Vec::new();
    for x in 1..=num_shares as i64 {
        let mut y = 0;
        let mut x_pow = 1;

        for &coeff in &coefficients {
            y = (y + coeff * x_pow) % prime;
            x_pow = (x_pow * x) % prime;
        }

        shares.push((x, y));
    }

    shares
}



fn reconstruct_secret(shares: &[(i64, i64)], prime: i64) -> i64 {
    let mut secret = 0;

    for (xj, yj) in shares {
        let mut num = 1;
        let mut den = 1;

        for (xm, _) in shares {
            if xm != xj {
                num = (num * (-xm + prime) % prime) % prime;
                den = (den * (xj - xm + prime) % prime) % prime;
            }
        }

        let lagrange_coeff = (num * mod_inv(den, prime)) % prime;
        secret = (secret + (*yj * lagrange_coeff) % prime) % prime;
    }

    return (secret + prime) % prime;
}

fn main() {
    let secret = 1234;
    let threshold = 5;
    let num_shares = 5;
    let prime = 1613;

    let shares = generate_shares(secret, threshold, num_shares, prime);
    println!("Shares: {:?}", shares);

    let subset = &shares[..threshold];
    let recovered_secret = reconstruct_secret(subset, prime);
    println!("Recovered Secret: {}", recovered_secret);
}
