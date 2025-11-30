# ADSIM Version Management & Release System

## ✅ Implementation Complete

This document summarizes the version management and release system implemented for ADSIM v0.1.0.

---

## 🎯 What Was Implemented

### 1. Centralized Version Management

**Files Created:**
- ✅ `VERSION` - Single source of truth (contains: `0.1.0`)
- ✅ `src/version.jl` - Julia module to read and export version
- ✅ `scripts/update_version.jl` - Automated version update script

**Key Functions:**
- `get_version()` - Returns version string from VERSION file
- `get_version_string()` - Returns formatted "ADSIM vX.Y.Z"

**Version Update Locations (9 files):**
1. `VERSION` - Central version file
2. `Project.toml` - Julia package version
3. `src/ADSIM.jl` - Header comment
4. `src/kernel.jl` - Header comment
5. `Problemtype/ADSIM_2025.gid/ADSIM_2025.tcl` - GiD splash screen
6. `Problemtype/ADSIM_2025.gid/ADSIM_2025.xml` - GiD metadata
7. `Problemtype/ADSIM_2025.gid/scripts/WriteMaterialData.tcl` - Output header
8. `Problemtype/ADSIM_2025.gid/scripts/WriteCalculationData.tcl` - Output header
9. `Problemtype/ADSIM_2025.gid/scripts/WriteMeshFile.tcl` - Output header

### 2. Version Display Integration

**Integrated in:**
- ✅ **Log file banner** (`src/kernel.jl` lines 99-103) - Shows version on every run
- ✅ **Command-line flag** (`julia kernel.jl --version`) - Quick version check
- ✅ **VTK output headers** (`src/write_vtk.jl`) - Version in ParaView files
- ✅ **GiD splash screen** - Version displayed on startup
- ✅ **GiD About dialog** - Version in interface

### 3. GitHub Actions Workflows

**Files Created:**
- ✅ `.github/workflows/release.yml` - Automated release builds on tags
- ✅ `.github/workflows/test.yml` - CI testing on push/PR

**Release Workflow Features:**
- Triggered automatically when pushing version tags (e.g., `v0.1.0`)
- Builds Windows executable using PackageCompiler.jl
- Validates version consistency (VERSION file vs git tag)
- Creates GitHub release with auto-generated notes
- Uploads executable as ZIP artifact
- Marks alpha/beta/rc versions as pre-releases

**Test Workflow Features:**
- Runs on every push to main/develop branches
- Tests multiple Julia versions (1.8, 1.9, 1.10)
- Tests on Ubuntu and Windows
- Runs example problems (Advection, Diffusion, Reaction)
- Validates version consistency across files
- Checks for unreplaced placeholders

### 4. Documentation

**Files Created:**
- ✅ `CHANGELOG.md` - Version history in Keep a Changelog format
- ✅ `.github/RELEASE_TEMPLATE.md` - Release checklist for maintainers
- ✅ `.github/HOW_TO_RELEASE.md` - Complete step-by-step release guide

**Files Updated:**
- ✅ `README.md` - Added version badge, installation instructions, release info

---

## 🚀 How to Use

### Checking Current Version

```bash
# From VERSION file
cat VERSION

# From command line
julia src/kernel.jl --version

# From Julia REPL
using ADSIM
println(get_version())
```

### Updating Version

```bash
# Update to new version (e.g., 0.2.0)
julia scripts/update_version.jl 0.2.0

# Verify changes
git diff

# All 9 files should be updated automatically
```

### Creating a Release

```bash
# 1. Update CHANGELOG.md with new version

# 2. Update version
julia scripts/update_version.jl 0.2.0

# 3. Commit changes
git add .
git commit -m "Bump version to v0.2.0"

# 4. Create annotated tag
git tag -a v0.2.0 -m "Release v0.2.0 - Brief description"

# 5. Push (triggers automated build!)
git push origin main
git push origin v0.2.0

# 6. Monitor workflow at:
# https://github.com/luisez1988/ADSIM/actions
```

---

## 📋 Version Update Verification

After running `julia scripts/update_version.jl 0.1.0`, verify these locations show `0.1.0`:

```bash
# Check all version locations
grep -r "0.1.0" VERSION
grep "version = " Project.toml
grep "# v" src/ADSIM.jl | head -1
grep "# v" src/kernel.jl | head -1
grep "adsim_version" Problemtype/ADSIM_2025.gid/ADSIM_2025.tcl
grep "Version" Problemtype/ADSIM_2025.gid/ADSIM_2025.xml
grep "ADSIM_version" Problemtype/ADSIM_2025.gid/scripts/*.tcl
```

---

## 🔍 File Structure Overview

```
ADSIM/
├── VERSION                                    # ← Central version file
├── Project.toml                               # ← Julia package metadata
├── CHANGELOG.md                               # ← Version history
├── README.md                                  # ← Updated with badges
│
├── src/
│   ├── version.jl                            # ← Version reader module
│   ├── ADSIM.jl                              # ← Exports version functions
│   ├── kernel.jl                             # ← Main with --version flag
│   └── write_vtk.jl                          # ← VTK with version header
│
├── scripts/
│   └── update_version.jl                     # ← Automated version updater
│
├── Problemtype/ADSIM_2025.gid/
│   ├── ADSIM_2025.tcl                        # ← Splash screen version
│   ├── ADSIM_2025.xml                        # ← Metadata version
│   └── scripts/
│       ├── WriteMaterialData.tcl             # ← Output file version
│       ├── WriteCalculationData.tcl          # ← Output file version
│       └── WriteMeshFile.tcl                 # ← Output file version
│
└── .github/
    ├── workflows/
    │   ├── release.yml                       # ← Automated release workflow
    │   └── test.yml                          # ← CI testing workflow
    ├── RELEASE_TEMPLATE.md                   # ← Release checklist
    └── HOW_TO_RELEASE.md                     # ← Complete release guide
```

---

## ✨ Features

### Automated Version Updates

The `update_version.jl` script automatically:
- ✅ Validates semantic versioning format
- ✅ Updates VERSION file
- ✅ Updates Project.toml
- ✅ Updates all Julia source headers
- ✅ Updates all GiD TCL/XML files
- ✅ Preserves file formatting
- ✅ Provides clear success feedback

### Version Display

Version appears in:
- ✅ Log file header: `Version: 0.1.0`
- ✅ VTK files: `ADSIM v0.1.0 - Time Step...`
- ✅ Command line: `julia kernel.jl --version` → `ADSIM v0.1.0`
- ✅ GiD splash screen: `ADSIM v0.1.0`
- ✅ Generated TOML/mesh files: `ADSIM_version = "2025 v0.1.0"`

### CI/CD Automation

**On every push to main/develop:**
- ✅ Tests run automatically
- ✅ Multiple Julia versions tested
- ✅ Example problems validated
- ✅ Version consistency checked

**On tag push (e.g., v0.1.0):**
- ✅ Executable builds automatically
- ✅ GitHub release created automatically
- ✅ ZIP artifact uploaded automatically
- ✅ Release notes extracted from CHANGELOG

---

## 🎓 Best Practices

### Semantic Versioning

Follow semver (MAJOR.MINOR.PATCH):
- **MAJOR** (1.0.0): Breaking changes, incompatible API
- **MINOR** (0.2.0): New features, backwards compatible
- **PATCH** (0.1.1): Bug fixes only

### Pre-release Labels

- **Alpha**: `0.1.0-alpha` - Early testing
- **Beta**: `0.2.0-beta` - Feature complete, testing
- **RC**: `1.0.0-rc.1` - Release candidate

### Workflow

1. Make changes
2. Update CHANGELOG.md
3. Run `julia scripts/update_version.jl X.Y.Z`
4. Commit and tag
5. Push tag
6. GitHub Actions handles the rest!

---

## 🔧 Troubleshooting

### Version Mismatch Error

**Problem:** VERSION file doesn't match git tag

**Solution:**
```bash
julia scripts/update_version.jl 0.1.0
git add VERSION
git commit --amend --no-edit
git tag -d v0.1.0
git tag -a v0.1.0 -m "Release v0.1.0"
git push --force origin main
git push origin v0.1.0
```

### Build Fails

**Problem:** PackageCompiler fails in GitHub Actions

**Solution:**
1. Test locally: `julia build_app.jl`
2. Check Project.toml and Manifest.toml are committed
3. Verify dependencies install: `julia -e 'using Pkg; Pkg.instantiate()'`
4. Check workflow logs for specific error

### Placeholder Found

**Problem:** CI detects "v0.x.x" placeholder

**Solution:**
```bash
# Find remaining placeholders
grep -r "v0\.x\.x" --include="*.jl" --include="*.tcl" --include="*.xml" .

# Re-run version update
julia scripts/update_version.jl 0.1.0
```

---

## 📊 Testing Checklist

Before creating v0.1.0 release:

- [x] VERSION file contains `0.1.0`
- [x] All 9 files updated to v0.1.0
- [x] `julia kernel.jl --version` returns `ADSIM v0.1.0`
- [x] Version appears in log file banner
- [ ] Build executable locally: `julia build_app.jl`
- [ ] Test executable: `.\ADSIM_app\bin\ADSIM.exe --version`
- [ ] Run test case: `julia -t8 kernel.jl Pani_etal`
- [ ] Verify VTK header shows v0.1.0
- [ ] CHANGELOG.md has v0.1.0 entry
- [ ] README.md updated with badges
- [ ] All tests pass: `julia --project=. -e 'using Pkg; Pkg.test()'`

---

## 🎉 Ready for First Release!

Everything is set up for your first GitHub release. Follow these steps:

1. **Build and test locally** (optional but recommended):
   ```bash
   julia build_app.jl
   .\ADSIM_app\bin\ADSIM.exe --version
   ```

2. **Commit all changes:**
   ```bash
   git add .
   git commit -m "Prepare for first release v0.1.0"
   ```

3. **Create tag:**
   ```bash
   git tag -a v0.1.0 -m "First alpha release of ADSIM"
   ```

4. **Push to GitHub:**
   ```bash
   git push origin main
   git push origin v0.1.0
   ```

5. **Monitor workflow:**
   - Visit: https://github.com/luisez1988/ADSIM/actions
   - Wait for green checkmark ✅

6. **Verify release:**
   - Visit: https://github.com/luisez1988/ADSIM/releases
   - Download and test ZIP file

---

## 📚 Additional Resources

- **Semantic Versioning:** https://semver.org/
- **Keep a Changelog:** https://keepachangelog.com/
- **GitHub Actions:** https://docs.github.com/en/actions
- **PackageCompiler.jl:** https://github.com/JuliaLang/PackageCompiler.jl

---

## 🤝 Contributing

When contributing version-related changes:

1. Never manually edit version strings
2. Always use `scripts/update_version.jl`
3. Update CHANGELOG.md for all notable changes
4. Test locally before pushing tags
5. Follow semantic versioning guidelines

---

**Implementation Date:** November 30, 2025  
**Current Version:** 0.1.0  
**Status:** ✅ Ready for Release

---

For questions or issues with the version management system, refer to `.github/HOW_TO_RELEASE.md` or open an issue on GitHub.
