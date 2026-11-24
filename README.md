# Basilisk-101: High-Fidelity Simulations Using Basilisk C

[![Issues](https://img.shields.io/github/issues/comphy-lab/Basilisk-101)](https://github.com/comphy-lab/Basilisk-101/issues)
[![License](https://img.shields.io/github/license/comphy-lab/Basilisk-101)](https://github.com/comphy-lab/Basilisk-101/blob/main/LICENSE)
[![Last Commit](https://img.shields.io/github/last-commit/comphy-lab/Basilisk-101)](https://github.com/comphy-lab/Basilisk-101/commits/main)
[![Course](https://img.shields.io/badge/Course-March%202025-blue)](https://comphy-lab.org/teaching/2025-Basilisk101-Madrid)
[![Basilisk](https://img.shields.io/badge/Basilisk-Compatible-green)](http://basilisk.fr/)
![CodeRabbit Pull Request Reviews](https://img.shields.io/coderabbit/prs/github/comphy-lab/Basilisk-101?utm_source=oss&utm_medium=github&utm_campaign=comphy-lab%2FBasilisk-101&labelColor=171717&color=FF570A&link=https%3A%2F%2Fcoderabbit.ai&label=CodeRabbit+Reviews)

Repository for the Basilisk C course at Universidad Carlos III de Madrid (March 10-13, 2025).

## Course Description

This intensive course provides a comprehensive introduction to high-fidelity simulations using Basilisk C, a powerful computational framework for fluid dynamics. Participants will learn to implement and solve complex fluid mechanics problems with an emphasis on multiphase flows, interface dynamics, and non-Newtonian rheology.

### What You'll Learn

- **Think before you compute!** Understanding the physics before implementation
- **Writing the first code in Basilisk C** Getting comfortable with the framework
- **Solving conservation equations** Numerical approaches to fluid dynamics
- **Interface tracking methods** Capturing multiphase phenomena accurately
- **Non-Newtonian flows** Modeling complex rheological behaviors

### Course Details

- **Dates**: March 10-13, 2025
- **Location**: Universidad Carlos III de Madrid, Spain
- **Duration**: 4 days, full-time
- **Prerequisites**:
  - Basic knowledge of fluid mechanics
  - Experience with programming (any language, C preferred)
  - Understanding of partial differential equations
  - Laptop with ability to compile C code

For complete course schedule and registration information, visit: [Course Website](https://comphy-lab.org/teaching/2025-Basilisk101-Madrid)

## Repository Structure

```
├-- basilisk/src/: Core Basilisk CFD library (reference only)
├── postProcess/: Project-specific post-processing tools
├── src-local/: Custom header files extending Basilisk functionality
├── testCases/: Test cases with their own Makefile
```

## Installation

Basilisk requires a C99-compliant compiler and GNU make. Installation can be done in two ways:

### Global Installation

#### Using darcs (recommended)
```bash
sudo apt install darcs make gawk
darcs clone https://basilisk.fr/basilisk
cd basilisk/src

# For Linux/Ubuntu users (preferred mode of operation)
ln -s config.gcc config

# For Mac users
# ln -s config.osx config

make
```

#### Using a tarball
```bash
wget http://basilisk.fr/basilisk/basilisk.tar.gz
tar xzf basilisk.tar.gz
cd basilisk/src

# For Linux/Ubuntu users (preferred mode of operation)
ln -s config.gcc config

# For Mac users
# ln -s config.osx config

make
```

#### Add to your shell configuration (.bashrc or .zshrc)
```bash
echo "export BASILISK=$PWD" >> ~/.bashrc
echo 'export PATH=$PATH:$BASILISK' >> ~/.bashrc
```
Or for zsh users:
```bash
echo "export BASILISK=$PWD" >> ~/.zshrc
echo 'export PATH=$PATH:$BASILISK' >> ~/.zshrc
```

### Repository Level Installation

For project-specific installations, you can use the provided `reset_install_requirements.sh` script which:
- Builds the Basilisk source committed in `basilisk/` (nothing is downloaded)
- Sets up environment variables locally (in `.project_config`)
- Automatically detects your OS (MacOS or Linux) and uses appropriate configuration
- Verifies the installation

#### Basic usage:
```bash
# Build the committed Basilisk source
./reset_install_requirements.sh

# Rebuild it from scratch (the source itself is not replaced)
./reset_install_requirements.sh --hard

# Load the environment settings for your current shell session
source .project_config
```

The script will create a `.project_config` file in your project root with the necessary environment variables. This approach avoids modifying your global shell configuration and keeps the Basilisk installation contained within your project.

### Windows Subsystem for Linux (WSL) Compatibility

Testing on WSL is currently incomplete. In principle, the Linux installation instructions should work for WSL environments. If you encounter any issues while installing or running Basilisk on WSL, please report them by [opening a bug report](https://github.com/comphy-lab/Basilisk-101/issues/new?template=bug_report.md&labels=bug,wsl).

### Complete Installation Instructions

For more detailed installation instructions, including configuration for different systems, setting up environment variables, installing additional dependencies, and optional libraries, please refer to the official installation guide at [http://basilisk.fr/src/INSTALL](http://basilisk.fr/src/INSTALL).

### Running the codes

To use the make file do:
```bash
CFLAGS=-DDISPLAY=-1 make NAME-of-File.tst
```

## Reporting Issues and Feature Requests

We use GitHub Issues to track bugs, feature requests, and example requests for this course. When creating an issue, please select the appropriate template to help us address your needs efficiently.

### Issue Templates

#### Bug Report:
[Report here](https://github.com/comphy-lab/Basilisk-101/issues/new?template=bug_report.md)

- For problems with installation, compilation, or running code. 
Please include:
- Detailed description of the issue
- Your environment (OS, compiler version)
- Steps to reproduce
- Expected vs. actual behavior
- Error messages and logs
- Code snippets or files that demonstrate the issue

#### Feature/Topic Request:
[Report here](https://github.com/comphy-lab/Basilisk-101/issues/new?template=feature_request.md)
- For requesting specific topics or functionality
- Coverage of specific topics in the course
- New examples or tutorials
- Additional functionality in the codebase
- Improvements to existing materials

#### Example Request:
[Report here](https://github.com/comphy-lab/Basilisk-101/issues/new?template=example_request.md)
- For requesting specific examples that demonstrate:
- Particular Basilisk features
- Solutions to common problems
- Implementation of specific physics or numerical methods

#### General Question:
[Report here](https://github.com/comphy-lab/Basilisk-101/issues/new?template=general_question.md)
- For any other questions

### How to Create an Issue

1. Go to the ["Issues" tab](https://github.com/comphy-lab/Basilisk-101/issues) in the GitHub repository
2. Click the ["New Issue"](https://github.com/comphy-lab/Basilisk-101/issues/new/choose) button
3. Select the appropriate template from the options
4. Fill in the required information according to the template
5. Add relevant labels if available
6. Submit the issue