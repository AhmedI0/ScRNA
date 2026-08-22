#!/usr/bin/env bash
# Step 09 - one .h5ad per library (Gene, GeneFull and spliced/unspliced layers).
source "$(dirname "$0")/common.sh"
for lib in $(library_ids); do
    solo="$SOLO_DIR/$lib/Solo.out"
    require_file "$solo/Gene/Summary.csv" "run 06_starsolo.sh"
    in_env "$ENV_UPSTREAM" python "$PKG_DIR/python/starsolo_to_h5ad.py" "$solo" "$lib" "$SOLO_DIR/h5ad/$lib.h5ad"
done
info "h5ad export finished: $SOLO_DIR/h5ad"
