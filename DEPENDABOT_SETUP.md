# Dependabot Configuration for Cramino

This document explains the Dependabot setup for automated dependency management in the cramino project.

## Files Added

### 1. `.github/dependabot.yml`

This is the main Dependabot configuration file that:

- Checks for Cargo dependency updates weekly on Mondays at 9:00 AM (Europe/Brussels timezone)
- Groups all patch updates together to reduce PR noise
- Assigns dependency PRs to @wdecoster for review
- Labels PRs with "dependencies" and "rust" tags
- Allows up to 10 open pull requests at once

### 2. `.github/workflows/dependabot.yml`

This workflow provides automated handling of Dependabot PRs:

- **Auto-merge**: Automatically merges patch and minor updates after tests pass
- **Manual review**: Flags major version updates for manual review
- **Safety checks**: Waits for the "build" job from the test workflow to succeed before merging

## How It Works

1. **Weekly Updates**: Every Monday, Dependabot scans the `Cargo.toml` file for outdated dependencies
2. **Pull Request Creation**: Creates PRs for available updates, grouping patch updates together
3. **Automated Testing**: The existing test workflow runs automatically on all Dependabot PRs
4. **Smart Merging**:
   - Patch/minor updates: Auto-merged after tests pass
   - Major updates: Flagged for manual review with a comment

## Configuration Details

- **Schedule**: Weekly on Mondays at 9:00 AM Europe/Brussels time
- **Timezone**: Configured for European working hours
- **Update Types**: Both direct and indirect dependencies
- **Grouping**: Patch updates are grouped to reduce noise
- **Review Process**: All PRs assigned to @wdecoster

## Benefits

- ✅ **Automated Security**: Quickly applies security patches
- ✅ **Reduced Maintenance**: Automatically keeps dependencies current
- ✅ **Safe Defaults**: Major updates require manual review
- ✅ **CI Integration**: Only merges when tests pass
- ✅ **Organized PRs**: Groups related updates together

## Monitoring

You can monitor Dependabot activity in the GitHub repository:

- **Security tab**: View security advisories and automated fixes
- **Pull requests**: See pending dependency updates
- **Actions tab**: Monitor the auto-merge workflow execution
