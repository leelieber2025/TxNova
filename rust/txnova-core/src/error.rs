use thiserror::Error;

#[derive(Debug, Error)]
pub enum CoreError {
    #[error("{0}")]
    Fail(String),
    #[error("{0}")]
    Io(#[from] std::io::Error),
    #[error("invalid JSON: {0}")]
    Json(#[from] serde_json::Error),
}

impl CoreError {
    pub fn fail(msg: impl Into<String>) -> Self {
        Self::Fail(msg.into())
    }
}

pub type Result<T> = std::result::Result<T, CoreError>;
