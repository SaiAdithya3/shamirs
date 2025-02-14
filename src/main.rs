use rand::Rng;
use std::error::Error;
use std::fmt;

#[derive(Debug)]
pub enum ShamirError {
    InvalidShares,
    InvalidParameters,
    ReconstructionError,
}

impl fmt::Display for ShamirError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            ShamirError::InvalidShares => write!(f, "Invalid number of shares"),
            ShamirError::InvalidParameters => write!(f, "Invalid parameters"),
            ShamirError::ReconstructionError => write!(f, "Failed to reconstruct secret"),
        }
    }
}

impl Error for ShamirError {}

pub struct Share {
    x: u64,
    y: u64,
}
#[derive(Debug)]
pub struct ShamirScheme {
    prime: u64,
    minimum_shares: usize,
    total_shares: usize,
}

impl ShamirScheme {
    pub fn new(minimum_shares: usize, total_shares: usize, prime: u64) -> Result<Self, ShamirError> {
        if minimum_shares > total_shares || minimum_shares == 0 {
            return Err(ShamirError::InvalidParameters);
        }
        
        Ok(ShamirScheme {
            prime,
            minimum_shares,
            total_shares,
        })
    }

    fn generate_polynomial(&self, secret: u64) -> Vec<u64> {
        let mut rng = rand::rng(); 
        let mut coefficients = vec![secret];
        
        for _ in 1..self.minimum_shares {
            coefficients.push(rng.random_range(1..self.prime));
        }
        
        coefficients
    }

    fn evaluate_polynomial(&self, coefficients: &[u64], x: u64) -> u64 {
        let mut result = 0;
        for &coeff in coefficients.iter().rev() {
            result = (result * x + coeff) % self.prime;
        }
        result
    }

    pub fn generate_shares(&self, secret: u64) -> Result<Vec<Share>, ShamirError> {
        if secret >= self.prime {
            return Err(ShamirError::InvalidParameters);
        }

        let polynomial = self.generate_polynomial(secret);
        let mut shares = Vec::with_capacity(self.total_shares);

        for x in 1..=self.total_shares {
            let y = self.evaluate_polynomial(&polynomial, x as u64);
            shares.push(Share {
                x: x as u64,
                y,
            });
        }

        Ok(shares)
    }

    fn mod_inverse(&self, mut a: i64, mut m: i64) -> Option<u64> {
        if m == 1 {
            return Some(0);
        }

        let m_orig = m;
        let mut x = 1i64;
        let mut y = 0i64;

        while a > 1 {
            let q = a / m;
            let t = m;
            m = a % m;
            a = t;
            let t = y;
            y = x - q * y;
            x = t;
        }

        if x < 0 {
            x += m_orig;
        }

        Some(x as u64)
    }

    pub fn reconstruct_secret(&self, shares: &[Share]) -> Result<u64, ShamirError> {
        if shares.len() < self.minimum_shares {
            return Err(ShamirError::InvalidShares);
        }

        let mut secret = 0u64;

        for i in 0..shares.len() {
            let mut numerator = 1i64;
            let mut denominator = 1i64;

            for j in 0..shares.len() {
                if i != j {
                    numerator = (numerator * -(shares[j].x as i64)) % self.prime as i64;
                    denominator = (denominator * ((shares[i].x as i64) - (shares[j].x as i64))) % self.prime as i64;
                }
            }

            if denominator < 0 {
                denominator += self.prime as i64;
            }

            let inverse = self.mod_inverse(denominator, self.prime as i64)
                .ok_or(ShamirError::ReconstructionError)?;

            let lagrange = (shares[i].y as u64 * numerator.unsigned_abs() as u64 * inverse) % self.prime;
            secret = (secret + lagrange) % self.prime;
        }

        Ok(secret)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_secret_sharing() {
        
        let prime = 2311u64;
        let scheme = ShamirScheme::new(3, 5, prime).unwrap();
        
        let secret = 1234u64;
        let shares = scheme.generate_shares(secret).unwrap();
        
        let min_shares = &shares[0..3];
        let reconstructed = scheme.reconstruct_secret(min_shares).unwrap();
        assert_eq!(secret, reconstructed);
        
        let reconstructed = scheme.reconstruct_secret(&shares).unwrap();
        assert_eq!(secret, reconstructed);
    }

    #[test]
    fn test_invalid_parameters() {
        let result = ShamirScheme::new(5, 3, 2311);
        assert!(result.is_err());
    }
}

fn main() {
    let scheme = ShamirScheme::new(8, 15, 2311).unwrap();
    
    let secret = 1234u64;
    
    match scheme.generate_shares(secret) {
        Ok(shares) => {
            println!("Generated shares: {:?}", shares.iter().map(|s| (s.x, s.y)).collect::<Vec<_>>());
            
             let min_shares = &shares[..8];
            match scheme.reconstruct_secret(min_shares) {
                Ok(recovered) => println!("Recovered secret: {}", recovered),
                Err(e) => println!("Failed to recover secret: {}", e),
            }
        },
        Err(e) => println!("Failed to generate shares: {}", e),
    }
    match scheme.generate_shares(secret) {
        Ok(shares) => {
            println!("Generated shares: {:?}", shares.iter().map(|s| (s.x, s.y)).collect::<Vec<_>>());
            
             let min_shares = &shares[..9];
            match scheme.reconstruct_secret(min_shares) {
                Ok(recovered) => println!("Recovered secret: {}", recovered),
                Err(e) => println!("Failed to recover secret: {}", e),
            }
        },
        Err(e) => println!("Failed to generate shares: {}", e),
    }
    match scheme.generate_shares(secret) {
        Ok(shares) => {
            println!("Generated shares: {:?}", shares.iter().map(|s| (s.x, s.y)).collect::<Vec<_>>());
            
             let min_shares = &shares[..10];
            match scheme.reconstruct_secret(min_shares) {
                Ok(recovered) => println!("Recovered secret: {}", recovered),
                Err(e) => println!("Failed to recover secret: {}", e),
            }
        },
        Err(e) => println!("Failed to generate shares: {}", e),
    }
}