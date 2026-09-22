#!/usr/bin/env bash
set -euo pipefail
if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "Usage: $0 OWNER/REPO [VERSION, default v0.0.37]" >&2
  echo "The repository must already exist and gh must be authenticated." >&2
  exit 2
fi
REPO="$1"
VERSION="${2:-v0.0.37}"
[[ "$REPO" =~ ^[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+$ ]] || exit 2
[[ "$VERSION" =~ ^v[0-9]+\.[0-9]+\.[0-9]+$ ]] || exit 2
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
gh auth status --hostname github.com
gh repo view "$REPO" --json nameWithOwner --jq .nameWithOwner
[[ -z "$(git status --porcelain)" ]] || { echo "Commit local changes first" >&2; exit 1; }
python -m unittest discover -s tests -v
if git rev-parse --verify "refs/tags/$VERSION" >/dev/null 2>&1; then
  [[ "$(git rev-parse "$VERSION^{commit}")" == "$(git rev-parse HEAD)" ]] || {
    echo "Tag does not match HEAD; refusing to replace a release tag" >&2; exit 1;
  }
else
  git tag "$VERSION"
fi
mkdir -p dist
git archive --format=tar.gz --prefix="profold-$VERSION/" \
  -o "dist/profold-$VERSION.tar.gz" "$VERSION"
(cd dist && sha256sum "profold-$VERSION.tar.gz" > "SHA256SUMS-$VERSION.txt")
# Authentication is supplied by gh; no token is embedded in the remote URL.
git -c credential.helper= -c 'credential.helper=!gh auth git-credential' \
  push "https://github.com/$REPO.git" HEAD:refs/heads/main "refs/tags/$VERSION"
gh release create "$VERSION" --repo "$REPO" --verify-tag \
  --title "$VERSION - ProFold" --notes-file RELEASE_NOTES.md \
  "dist/profold-$VERSION.tar.gz" "dist/SHA256SUMS-$VERSION.txt"
