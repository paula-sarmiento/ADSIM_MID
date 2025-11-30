# How to Create a GitHub Release for ADSIM

This guide walks you through creating your first GitHub release and all subsequent releases.

## Prerequisites

Before creating a release, ensure you have:

- ✅ Git installed and configured
- ✅ GitHub account with push access to the repository
- ✅ Julia 1.8+ installed (for local testing)
- ✅ All changes committed and tests passing

---

## Step-by-Step Guide

### Step 1: Update the CHANGELOG

Before releasing, document all changes in `CHANGELOG.md`:

```markdown
## [0.2.0] - 2025-12-15

### Added
- New feature X
- Support for Y

### Fixed
- Bug in Z function
- Memory leak in solver

### Changed
- Improved performance of W module
```

### Step 2: Update the Version

Use the automated version update script:

```bash
# Navigate to project root
cd c:\Users\zamcr\Dcuments\GitHub\ADSIM

# Update version (example: 0.2.0)
julia scripts/update_version.jl 0.2.0
```

This automatically updates:
- ✅ VERSION file
- ✅ Project.toml
- ✅ All Julia source headers
- ✅ GiD TCL and XML files
- ✅ TCL script files

**Verify the changes:**
```bash
git diff
```

You should see version updates in 9 files.

### Step 3: Test Locally

**Build and test the executable:**

```bash
# Build standalone executable
julia build_app.jl

# Test the version
.\ADSIM_app\bin\ADSIM.exe --version

# Run a test case
cd src
julia -t8 kernel.jl Pani_etal
```

**Verify output:**
- Check that version appears in log file: `src/output/Pani_etal.log`
- Check VTK file headers contain correct version
- Ensure no errors during execution

### Step 4: Commit and Tag

**Commit the version changes:**

```bash
git add .
git commit -m "Bump version to v0.2.0"
```

**Create an annotated tag:**

```bash
# Format: git tag -a v<VERSION> -m "<MESSAGE>"
git tag -a v0.2.0 -m "Release v0.2.0 - Added feature X, fixed bug Y"
```

**View your tags:**
```bash
git tag -l
```

### Step 5: Push to GitHub

**Push the commit and tag:**

```bash
# Push the main branch
git push origin main

# Push the tag (this triggers the release workflow!)
git push origin v0.2.0
```

### Step 6: Monitor the GitHub Actions Workflow

1. **Go to GitHub Actions:**
   - Navigate to: `https://github.com/luisez1988/ADSIM/actions`
   - Or click the "Actions" tab in your repository

2. **Watch the "Release Build" workflow:**
   - You should see a new workflow run for your tag
   - It will show real-time progress through these steps:
     - ✅ Checkout code
     - ✅ Setup Julia
     - ✅ Install dependencies
     - ✅ Build executable with PackageCompiler
     - ✅ Package as ZIP
     - ✅ Create GitHub release
     - ✅ Upload artifacts

3. **Wait for completion:**
   - Typically takes 10-15 minutes for the build
   - A green checkmark ✅ means success
   - A red X ❌ means failure (check logs)

### Step 7: Verify the Release

1. **Go to Releases page:**
   - Navigate to: `https://github.com/luisez1988/ADSIM/releases`
   - Or click "Releases" on the right side of the repository

2. **Check your new release:**
   - Release title should be `v0.2.0`
   - Release notes should be auto-populated from CHANGELOG.md
   - Download the ZIP file: `ADSIM-v0.2.0-windows-x64.zip`

3. **Test the released executable:**
   ```bash
   # Extract ZIP to a test location
   # Run the executable
   .\ADSIM_app\bin\ADSIM.exe --version
   ```

---

## What Happens Automatically?

When you push a tag to GitHub, the workflow automatically:

1. ✅ **Validates version consistency** - Checks that VERSION file matches the tag
2. ✅ **Sets up Julia environment** - Installs Julia 1.8
3. ✅ **Installs dependencies** - Runs `Pkg.instantiate()`
4. ✅ **Builds executable** - Uses PackageCompiler to create `ADSIM.exe`
5. ✅ **Packages release** - Creates ZIP file with all necessary files
6. ✅ **Creates GitHub release** - Publishes release with auto-generated notes
7. ✅ **Uploads artifacts** - Attaches ZIP file to the release
8. ✅ **Marks as pre-release** - If version contains "alpha", "beta", or "rc"

---

## Release Naming Convention

Follow semantic versioning:

| Version Type | Format | Example | When to Use |
|--------------|--------|---------|-------------|
| **Major** | X.0.0 | v1.0.0 | Breaking changes, incompatible API |
| **Minor** | 0.X.0 | v0.2.0 | New features, backwards compatible |
| **Patch** | 0.0.X | v0.1.1 | Bug fixes only |
| **Alpha** | X.Y.Z-alpha | v0.2.0-alpha | Early testing |
| **Beta** | X.Y.Z-beta | v0.2.0-beta | Feature complete, testing |
| **RC** | X.Y.Z-rc.N | v1.0.0-rc.1 | Release candidate |

---

## Troubleshooting

### Problem: Workflow fails with "Version mismatch"

**Solution:**
```bash
# Ensure VERSION file matches your tag
cat VERSION  # Should show: 0.2.0 (without 'v')

# If mismatch, fix and re-tag:
julia scripts/update_version.jl 0.2.0
git add VERSION
git commit --amend --no-edit
git tag -d v0.2.0  # Delete old tag
git tag -a v0.2.0 -m "Release v0.2.0"
git push origin main --force
git push origin v0.2.0
```

### Problem: Build fails in GitHub Actions

**Solution:**
1. Check workflow logs for specific error
2. Test build locally: `julia build_app.jl`
3. Verify `Project.toml` and `Manifest.toml` are committed
4. Check that all dependencies install correctly

### Problem: Release created but executable is missing

**Solution:**
1. Check that build step completed successfully
2. Look for "Upload executable artifact" step in logs
3. Verify ZIP file was created correctly
4. Re-run the workflow from GitHub Actions UI

### Problem: Need to delete/fix a release

**Delete the release:**
```bash
# Delete tag locally
git tag -d v0.2.0

# Delete tag on GitHub
git push origin :refs/tags/v0.2.0
```

**Delete release on GitHub:**
1. Go to Releases page
2. Click on the release
3. Click "Delete" button
4. Then create a new corrected release

---

## First-Time Release (v0.1.0)

For your first release:

1. **All files are already at v0.1.0** (we just ran the update script)

2. **Commit everything:**
   ```bash
   git add .
   git commit -m "Prepare for first release v0.1.0
   
   - Add VERSION file and version management system
   - Add automated version update script
   - Add GitHub Actions workflows for CI/CD
   - Add CHANGELOG.md
   - Update README with installation instructions
   - Integrate version display in application"
   ```

3. **Create the tag:**
   ```bash
   git tag -a v0.1.0 -m "First alpha release of ADSIM
   
   This release includes:
   - Core finite element solver for multi-gas transport
   - GiD problem type for pre-processing
   - VTK output for ParaView visualization
   - Example test cases
   - Windows executable
   
   This is an alpha release for testing and feedback."
   ```

4. **Push to GitHub:**
   ```bash
   git push origin main
   git push origin v0.1.0
   ```

5. **Wait for workflow and verify release!**

---

## Manual Release (If Workflow Fails)

If GitHub Actions doesn't work, you can create a release manually:

1. **Build locally:**
   ```bash
   julia build_app.jl
   ```

2. **Create ZIP:**
   ```powershell
   $version = Get-Content VERSION
   Compress-Archive -Path ADSIM_app/* -DestinationPath "ADSIM-v$version-windows-x64.zip"
   ```

3. **Create release on GitHub:**
   - Go to: `https://github.com/luisez1988/ADSIM/releases/new`
   - Choose your tag: `v0.1.0`
   - Release title: `v0.1.0`
   - Copy release notes from CHANGELOG.md
   - Attach ZIP file
   - Check "This is a pre-release" for alpha/beta
   - Click "Publish release"

---

## After Release

1. ✅ Announce on relevant channels (if applicable)
2. ✅ Update documentation website (if applicable)
3. ✅ Start working on next version
4. ✅ Update CHANGELOG with "Unreleased" section

**Example for next version:**
```markdown
## [Unreleased]

### Added
- Work in progress...

### Fixed
- Fixes pending...
```

---

## Quick Reference Commands

```bash
# Check current version
cat VERSION

# Update version
julia scripts/update_version.jl X.Y.Z

# Build and test locally
julia build_app.jl
.\ADSIM_app\bin\ADSIM.exe --version

# Create release
git add .
git commit -m "Bump version to vX.Y.Z"
git tag -a vX.Y.Z -m "Release vX.Y.Z - Brief description"
git push origin main
git push origin vX.Y.Z

# Check status
# Visit: https://github.com/luisez1988/ADSIM/actions
```

---

## Getting Help

If you encounter issues:

1. **Check workflow logs** on GitHub Actions
2. **Test locally** to isolate the problem
3. **Review this guide** for common solutions
4. **Check GitHub Actions documentation** for workflow issues
5. **Open an issue** if you find a bug in the release process

---

**You're now ready to create your first release! 🎉**

Good luck with v0.1.0!
