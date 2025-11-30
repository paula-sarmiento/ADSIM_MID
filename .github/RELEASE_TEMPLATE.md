## Release Checklist for ADSIM

### Pre-Release Steps

- [ ] All tests passing locally
- [ ] All example problems run successfully
- [ ] Documentation updated for new features
- [ ] CHANGELOG.md updated with all changes
- [ ] Version number decided (follow semantic versioning)

### Version Update

1. **Update version using the automated script:**
   ```bash
   julia scripts/update_version.jl X.Y.Z
   ```
   
   This updates all 9 locations:
   - VERSION file
   - Project.toml
   - src/ADSIM.jl
   - src/kernel.jl
   - Problemtype/ADSIM_2025.gid/ADSIM_2025.tcl
   - Problemtype/ADSIM_2025.gid/ADSIM_2025.xml
   - Problemtype/ADSIM_2025.gid/scripts/*.tcl (3 files)

2. **Verify changes:**
   ```bash
   git diff
   ```

### Build and Test

- [ ] Build executable locally:
  ```bash
  julia build_app.jl
  ```

- [ ] Test executable:
  ```bash
  cd ADSIM_app/bin
  ./ADSIM.exe --version
  ./ADSIM.exe Pani_etal
  ```

- [ ] Verify output files generated correctly

### Commit and Tag

1. **Commit version changes:**
   ```bash
   git add .
   git commit -m "Bump version to vX.Y.Z"
   ```

2. **Create annotated tag:**
   ```bash
   git tag -a vX.Y.Z -m "Release vX.Y.Z - Brief description"
   ```

3. **Push to GitHub:**
   ```bash
   git push origin main
   git push origin vX.Y.Z
   ```

### GitHub Release

The GitHub Actions workflow will automatically:
- Build the Windows executable
- Create a GitHub release
- Upload the executable as a ZIP file
- Extract release notes from CHANGELOG.md

**Monitor the workflow:**
- Go to: https://github.com/luisez1988/ADSIM/actions
- Check that the "Release Build" workflow completes successfully

### Post-Release

- [ ] Verify release appears on: https://github.com/luisez1988/ADSIM/releases
- [ ] Download and test the released executable
- [ ] Update README badges if needed
- [ ] Announce release (if applicable)
- [ ] Update documentation website (if applicable)

### Semantic Versioning Guide

- **MAJOR** (X.0.0): Incompatible API changes
- **MINOR** (0.X.0): New functionality, backwards compatible
- **PATCH** (0.0.X): Bug fixes, backwards compatible

### Pre-release Labels

- **Alpha**: `vX.Y.Z-alpha` or `vX.Y.Z-alpha.N`
- **Beta**: `vX.Y.Z-beta` or `vX.Y.Z-beta.N`
- **Release Candidate**: `vX.Y.Z-rc` or `vX.Y.Z-rc.N`

### Rollback Procedure

If a release needs to be rolled back:

1. **Delete the tag locally and remotely:**
   ```bash
   git tag -d vX.Y.Z
   git push origin :refs/tags/vX.Y.Z
   ```

2. **Delete the GitHub release:**
   - Go to releases page
   - Click on the release
   - Click "Delete this release"

3. **Revert version changes:**
   ```bash
   git revert <commit-hash>
   git push origin main
   ```

### Troubleshooting

**Build fails in GitHub Actions:**
- Check the workflow logs
- Verify Project.toml and Manifest.toml are up to date
- Test the build locally first

**Version mismatch error:**
- Ensure VERSION file matches the git tag
- Re-run `julia scripts/update_version.jl X.Y.Z`

**Executable doesn't run:**
- Check Julia version compatibility
- Verify all dependencies are included in build
- Test on a clean Windows machine
