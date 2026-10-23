import math
import matplotlib.pyplot as plt

# ==============================================================================
# PHYSICAL CONSTANTS (CODATA values in SI units)
# ==============================================================================
h = 6.62607015e-34    # Planck's constant (J*s)
c = 2.99792458e8       # Speed of light in vacuum (m/s)
ev = 1.602176634e-19   # Elementary charge / 1 electronvolt in Joules (J)

# Sensible physical limits (keeps calculations within safe floating-point range)
MIN_LENGTH, MAX_LENGTH = 1e-18, 1e12      # metres
MIN_MASS, MAX_MASS = 1e-40, 1e10          # kilograms
MAX_N = 10**6                             # largest supported quantum number

# FIX: largest quantum numbers for which the graphs are still reliable.
# Above these, the graph would be aliased (too few points per lobe) or would
# need millions of nodal lines, so the graph is skipped. Numbers are unaffected.
PLOT_MAX_N_1D = 2500
PLOT_MAX_N_2D = 30
PLOT_MAX_N_3D = 6


# ==============================================================================
# SAFE INPUT HELPERS
# ==============================================================================

def ask_float(prompt):
    """Safely prompt the user for a finite floating-point number."""
    while True:
        try:
            x = float(input(prompt))
        except ValueError:
            print("Please enter a valid number.")
            continue
        if math.isfinite(x):
            return x
        print("Value must be finite.")


def ask_int(prompt):
    """Safely prompt the user for a positive integer quantum number."""
    while True:
        try:
            n = int(input(prompt))
        except ValueError:
            print("Please enter a whole number (1, 2, 3, ...).")
            continue
        if 1 <= n <= MAX_N:
            return n
        print("Quantum number must be between 1 and", MAX_N)


# ==============================================================================
# LENGTH INPUT AND UNIT CONVERSION
# ==============================================================================

def length(message, allow_zero=False):
    """
    Prompt for a physical length with standard unit conversion to metres.
    Supports Angstroms, Picometres, Nanometres, Micrometres, Millimetres,
    Centimetres, and Metres.
    """
    units = {
        "1": 1e-10,    # Angstrom (A)
        "2": 1e-12,    # Picometre (pm)
        "3": 1e-9,     # Nanometre (nm)
        "4": 1e-6,     # Micrometre (um)
        "5": 1e-3,     # Millimetre (mm)
        "6": 1e-2,     # Centimetre (cm)
        "7": 1.0       # Metre (m)
    }

    while True:
        print("\nEnter", message)
        x = ask_float("Enter numerical value: ")

        print("\nSelect unit:")
        print("1. Angstrom (A)")
        print("2. Picometre (pm)")
        print("3. Nanometre (nm)")
        print("4. Micrometre (um)")
        print("5. Millimetre (mm)")
        print("6. Centimetre (cm)")
        print("7. Metre (m)")

        u = input("Enter choice (1-7): ").strip()

        if u not in units:
            print("Invalid unit selection. Please choose 1 to 7.")
            continue

        x = x * units[u]

        if (allow_zero and x == 0) or MIN_LENGTH <= x <= MAX_LENGTH:
            print("Value in metres =", x)   # FIX: wording (works for positions too)
            return x

        print("Length must be between", MIN_LENGTH, "m and", MAX_LENGTH, "m.")


# ==============================================================================
# PARTICLE MASS INPUT AND UNIT CONVERSION
# ==============================================================================

def mass():
    """
    Prompt the user to select a known subatomic particle or specify
    a custom mass with standard unit conversions to kilograms.
    """
    particles = {
        "1": ("Electron", 9.1093837e-31),
        "2": ("Proton", 1.6726219e-27),
        "3": ("Neutron", 1.6749275e-27),
        "4": ("Alpha particle", 6.6446573e-27)
    }

    units = {
        "1": 1e-9,     # Microgram (ug) -> kg
        "2": 1e-6,     # Milligram (mg) -> kg
        "3": 1e-3,     # Gram (g) -> kg
        "4": 1.0       # Kilogram (kg)
    }

    while True:
        print("\n===================================")
        print("       PARTICLE MASS INPUT")
        print("===================================")
        print("1. Electron")
        print("2. Proton")
        print("3. Neutron")
        print("4. Alpha particle")
        print("5. Custom mass")

        choice = input("Enter your choice (1-5): ").strip()

        if choice in particles:
            name, m = particles[choice]
            print("Selected particle:", name)
            print("Mass =", format(m, ".6e"), "kg")
            return m

        if choice == "5":
            print("\nEnter the custom particle mass.")
            m = ask_float("Enter numerical value: ")

            print("\nSelect unit:")
            print("1. Microgram (ug)")
            print("2. Milligram (mg)")
            print("3. Gram (g)")
            print("4. Kilogram (kg)")

            u = input("Enter choice (1-4): ").strip()

            if u not in units:
                print("Invalid mass unit. Please choose 1 to 4.")
                continue

            m = m * units[u]

            if MIN_MASS <= m <= MAX_MASS:
                print("Custom mass =", format(m, ".6e"), "kg")
                return m

            print("Mass must be between", MIN_MASS, "kg and", MAX_MASS, "kg.")
            continue

        print("Invalid choice. Please select 1 to 5.")


# ==============================================================================
# ENERGY CALCULATIONS (1D, 2D, and 3D)
# ==============================================================================

def energy(state, m, L):
    """
    Calculate the total energy of a particle in an infinite potential well.
    Formula:
      E = (h^2 / (8 * m)) * sum( (n_i / L_i)^2 )
    Works for 1D, 2D, and 3D systems where 'state' is a tuple of quantum numbers
    and 'L' is a tuple of well lengths along each dimension.
    """
    total = 0.0
    for i in range(len(state)):
        total = total + (state[i]**2) / (L[i]**2)

    E = (h**2 * total) / (8.0 * m)
    return E


def energy_1d(n, m, L):
    """
    1D Infinite Potential Well Energy:
      E_n = (n^2 * h^2) / (8 * m * L^2)
    """
    return energy((n,), m, (L,))


def energy_2d(nx, ny, m, Lx, Ly):
    """
    2D Infinite Rectangular Potential Well Energy:
      E_nx,ny = (h^2 / (8 * m)) * ( (nx / Lx)^2 + (ny / Ly)^2 )
    """
    return energy((nx, ny), m, (Lx, Ly))


def energy_difference(initial, final, m, L):
    """
    Calculate the signed energy difference: delta_E = E_final - E_initial.
    Sums (n_final^2 - n_initial^2) / L^2 first to prevent catastrophic
    cancellation when subtracting very close or degenerate energy states.
    """
    terms = []
    for i in range(len(initial)):
        terms.append((final[i]**2 - initial[i]**2) / (L[i]**2))

    total = math.fsum(terms)
    size = math.fsum(abs(t) for t in terms)

    # If the result is negligible compared to the terms, it is numerical noise
    if abs(total) < 1e-14 * size:
        return 0.0

    return (h**2 * total) / (8.0 * m)


# ==============================================================================
# WAVEFUNCTION VALUES AT A POINT (1D, 2D, 3D)
# ==============================================================================

def psi_value(x, L, n, A):
    """
    Calculate the 1D wavefunction value at position x:
      psi(x) = A * sin(n * pi * x / L)
    Enforces boundary condition psi(0) = psi(L) = 0 and exact zeros at nodes.
    """
    if x <= 0 or x >= L or math.isclose(x, L, rel_tol=1e-12):
        return 0.0

    t = n * x / L

    # Nodes (where t is an integer) are exactly zero
    if abs(t - round(t)) < 1e-12:
        return 0.0

    return A * math.sin(math.pi * t)


def psi_value_2d(x, y, length_x, length_y, quantum_number_x, quantum_number_y, A):
    """
    Calculate the 2D wavefunction value at position (x, y):
      psi(x, y) = A * sin(nx * pi * x / Lx) * sin(ny * pi * y / Ly)
    where A = sqrt(4 / (Lx * Ly)).
    Vanishes at well boundaries and along internal nodal lines.
    """
    # Boundary conditions: wavefunction is zero at the well walls
    if x <= 0 or x >= length_x or math.isclose(x, length_x, rel_tol=1e-12):
        return 0.0
    if y <= 0 or y >= length_y or math.isclose(y, length_y, rel_tol=1e-12):
        return 0.0

    tx = quantum_number_x * x / length_x
    ty = quantum_number_y * y / length_y

    # Internal nodal lines are exactly zero
    if abs(tx - round(tx)) < 1e-12 or abs(ty - round(ty)) < 1e-12:
        return 0.0

    return A * math.sin(math.pi * tx) * math.sin(math.pi * ty)


def psi_value_3d(x, y, z, length_x, length_y, length_z, quantum_number_x, quantum_number_y, quantum_number_z, A):
    """
    Calculate the 3D wavefunction value at position (x, y, z):
      psi(x, y, z) = A * sin(nx * pi * x / Lx) * sin(ny * pi * y / Ly) * sin(nz * pi * z / Lz)
    where A = sqrt(8 / (Lx * Ly * Lz)).
    Vanishes at well boundaries and across internal nodal planes.
    """
    # Boundary conditions: wavefunction is zero at the well walls
    if x <= 0 or x >= length_x or math.isclose(x, length_x, rel_tol=1e-12):
        return 0.0
    if y <= 0 or y >= length_y or math.isclose(y, length_y, rel_tol=1e-12):
        return 0.0
    if z <= 0 or z >= length_z or math.isclose(z, length_z, rel_tol=1e-12):
        return 0.0

    tx = quantum_number_x * x / length_x
    ty = quantum_number_y * y / length_y
    tz = quantum_number_z * z / length_z

    # Internal nodal planes are exactly zero
    if abs(tx - round(tx)) < 1e-12 or abs(ty - round(ty)) < 1e-12 or abs(tz - round(tz)) < 1e-12:
        return 0.0

    return A * math.sin(math.pi * tx) * math.sin(math.pi * ty) * math.sin(math.pi * tz)


# ==============================================================================
# PROBABILITY CALCULATIONS IN A SUB-REGION
# ==============================================================================

def probability_between(n, a, b):
    """
    Analytical integration of probability density in 1D between a*L and b*L:
      integral_{a*L}^{b*L} (2/L) * sin^2(n * pi * x / L) dx
      = (b - a) - [sin(2*n*pi*b) - sin(2*n*pi*a)] / (2 * n * pi)
    Clamped to [0.0, 1.0] to eliminate floating-point noise.
    """
    p = (b - a) - (
        math.sin(2.0 * n * math.pi * b) -
        math.sin(2.0 * n * math.pi * a)
    ) / (2.0 * n * math.pi)

    return min(1.0, max(0.0, p))


def probability_between_2d(nx, ny, ax, bx, ay, by):
    """
    Calculate the probability in a 2D rectangular sub-region
    [ax*Lx, bx*Lx] x [ay*Ly, by*Ly] using separable 1D integrals:
      P = P_x * P_y
    """
    px = probability_between(nx, ax, bx)
    py = probability_between(ny, ay, by)
    return px * py


def probability_between_3d(nx, ny, nz, ax, bx, ay, by, az, bz):
    """
    Calculate the probability in a 3D rectangular sub-volume
    [ax*Lx, bx*Lx] x [ay*Ly, by*Ly] x [az*Lz, bz*Lz] using separable 1D integrals:
      P = P_x * P_y * P_z
    """
    px = probability_between(nx, ax, bx)
    py = probability_between(ny, ay, by)
    pz = probability_between(nz, az, bz)
    return px * py * pz


# ==============================================================================
# ENERGY DISPLAY IN JOULES AND ELECTRONVOLTS
# ==============================================================================

def show_energy(E):
    """Print an energy value formatted in Joules (J) and electronvolts (eV)."""
    print("Energy =", format(E, ".6e"), "J")
    print("Energy =", format(E / ev, ".6f"), "eV")


# ==============================================================================
# 1D WAVEFUNCTION ANALYSIS
# ==============================================================================

def wavefunction_1d():
    """Perform complete 1D potential well wavefunction analysis and plotting."""
    print("\n===================================")
    print("     1D WAVEFUNCTION ANALYSIS")
    print("===================================")
    print("Given wavefunction:")
    print("psi(x) = A sin(n*pi*x/L)")

    m = mass()
    L = length("the width of the potential well L")

    n = ask_int("\nEnter quantum number n (1, 2, 3, ...): ")

    # Normalization constant: A = sqrt(2 / L)
    A = math.sqrt(2.0 / L)
    print("\nNormalization constant A =", format(A, ".6e"))

    # Energy
    print("\n---------- PARTICLE ENERGY ----------")
    E = energy((n,), m, (L,))
    show_energy(E)

    # Probability density at a point
    print("\n---------- PROBABILITY DENSITY ----------")
    while True:
        x = length("the position x (0 <= x <= L)", allow_zero=True)
        if x <= L or math.isclose(x, L, rel_tol=1e-12):
            break
        print("Position must be between 0 and L.")

    psi = psi_value(x, L, n, A)
    density = psi**2

    print("Wavefunction psi(x) =", format(psi, ".6e"))
    print("Probability density =", format(density, ".6e"), "per metre")

    # Probability in a region
    print("\n---------- PROBABILITY IN A REGION ----------")
    print("Enter limits as fractions of L.")
    print("Example: 0.2 means position 0.2L.")

    while True:
        a = ask_float("Enter lower limit a/L (0 to 1): ")
        b = ask_float("Enter upper limit b/L (0 to 1): ")

        if 0 <= a < b <= 1:
            break
        print("Invalid limits. Require 0 <= a/L < b/L <= 1.")

    probability = probability_between(n, a, b)
    print("Probability =", format(probability, ".6f"))
    print("Probability percentage =", format(probability * 100.0, ".4f"), "%")

    # Probability density graph
    print("\n---------- PROBABILITY DENSITY GRAPH ----------")

    # FIX: skip the graph when n is too large to draw reliably
    if n > PLOT_MAX_N_1D:
        print("n is too large to draw a reliable graph (limit n <=", PLOT_MAX_N_1D, ").")
        print("Graph skipped. The numerical results above are still valid.")
        return

    print("Generating graph...")

    # 40 points per lobe is plenty (max 100000 points at n = 2500)
    points = max(500, 40 * n)

    x_values = []
    density_values = []

    for i in range(points):
        x = L * i / (points - 1)
        psi_pt = psi_value(x, L, n, A)
        x_values.append(x)
        density_values.append(psi_pt**2)

    plt.figure()
    plt.plot(x_values, density_values)
    plt.xlabel("Position x (m)")
    plt.ylabel("Probability density |psi(x)|^2 (m^-1)")
    plt.title("Probability Density for n = " + str(n))
    plt.grid(True)
    plt.show()


# ==============================================================================
# 2D WAVEFUNCTION ANALYSIS
# ==============================================================================

def wavefunction_2d():
    """Perform complete 2D potential well wavefunction analysis and 2D plotting."""
    print("\n===================================")
    print("     2D WAVEFUNCTION ANALYSIS")
    print("===================================")
    print("Given 2D wavefunction:")
    print("psi(x, y) = A * sin(nx*pi*x/Lx) * sin(ny*pi*y/Ly)")
    print("where A = sqrt(4 / (Lx * Ly))")

    m = mass()
    length_x = length("the well length Lx")
    length_y = length("the well length Ly")

    print("\nEnter quantum numbers (nx, ny):")
    quantum_number_x = ask_int("Enter quantum number nx (1, 2, 3, ...): ")
    quantum_number_y = ask_int("Enter quantum number ny (1, 2, 3, ...): ")

    # Normalization constant: A = sqrt(4 / (Lx * Ly))
    A = math.sqrt(4.0 / (length_x * length_y))
    print("\nNormalization constant A =", format(A, ".6e"))

    # Energy calculation: E = (h^2 / (8*m)) * (nx^2/Lx^2 + ny^2/Ly^2)
    print("\n---------- PARTICLE ENERGY ----------")
    energy_val = energy((quantum_number_x, quantum_number_y), m, (length_x, length_y))
    show_energy(energy_val)

    # Probability density at a point (x, y)
    print("\n---------- PROBABILITY DENSITY AT A POINT ----------")
    while True:
        x = length("the position x (0 <= x <= Lx)", allow_zero=True)
        if x <= length_x or math.isclose(x, length_x, rel_tol=1e-12):
            break
        print("Position x must be between 0 and Lx.")

    while True:
        y = length("the position y (0 <= y <= Ly)", allow_zero=True)
        if y <= length_y or math.isclose(y, length_y, rel_tol=1e-12):
            break
        print("Position y must be between 0 and Ly.")

    psi = psi_value_2d(x, y, length_x, length_y, quantum_number_x, quantum_number_y, A)
    probability_density = psi**2

    print("Wavefunction psi(x, y) =", format(psi, ".6e"))
    print("Probability density =", format(probability_density, ".6e"), "per m^2")

    # Probability in a region [ax*Lx, bx*Lx] x [ay*Ly, by*Ly]
    print("\n---------- PROBABILITY IN A REGION ----------")
    print("Enter limits along x as fractions of Lx (0 to 1).")
    print("Example: 0.2 means position 0.2 Lx.")
    while True:
        ax = ask_float("Enter lower limit ax/Lx (0 to 1): ")
        bx = ask_float("Enter upper limit bx/Lx (0 to 1): ")
        if 0 <= ax < bx <= 1:
            break
        print("Invalid limits. Require 0 <= ax/Lx < bx/Lx <= 1.")

    print("\nEnter limits along y as fractions of Ly (0 to 1).")
    print("Example: 0.2 means position 0.2 Ly.")
    while True:
        ay = ask_float("Enter lower limit ay/Ly (0 to 1): ")
        by = ask_float("Enter upper limit by/Ly (0 to 1): ")
        if 0 <= ay < by <= 1:
            break
        print("Invalid limits. Require 0 <= ay/Ly < by/Ly <= 1.")

    total_prob = probability_between_2d(quantum_number_x, quantum_number_y, ax, bx, ay, by)
    print("Probability =", format(total_prob, ".6f"))
    print("Probability percentage =", format(total_prob * 100.0, ".4f"), "%")

    # 2D Probability density graph
    print("\n---------- 2D PROBABILITY DENSITY GRAPH ----------")

    # FIX: skip the graph when nx or ny is too large to draw reliably
    if quantum_number_x > PLOT_MAX_N_2D or quantum_number_y > PLOT_MAX_N_2D:
        print("Quantum numbers are too large to draw a reliable graph (limit n <=", PLOT_MAX_N_2D, ").")
        print("Graph skipped. The numerical results above are still valid.")
        return

    print("Generating 2D graph...")

    # Grid resolution adapted to quantum numbers
    points_x = min(max(100, 30 * quantum_number_x), 300)
    points_y = min(max(100, 30 * quantum_number_y), 300)

    density_grid = []
    for j in range(points_y):
        y_val = length_y * j / (points_y - 1)
        row = []
        for i in range(points_x):
            x_val = length_x * i / (points_x - 1)
            psi_pt = psi_value_2d(
                x_val, y_val, length_x, length_y, quantum_number_x, quantum_number_y, A
            )
            row.append(psi_pt**2)
        density_grid.append(row)

    fig, ax = plt.subplots(figsize=(8, 6))
    im = ax.imshow(
        density_grid,
        extent=[0, length_x, 0, length_y],
        origin="lower",
        cmap="plasma",
        aspect="auto",
        vmin=0.0,
        vmax=A**2
    )
    cb = fig.colorbar(im, ax=ax)
    cb.set_label("Probability density |psi(x, y)|^2 (m^-2)")

    # Delineate internal nodal lines where probability density is zero (rho = 0)
    has_2d_nodes = False
    if quantum_number_x > 1:
        has_2d_nodes = True
        for kn in range(1, quantum_number_x):
            ax.axvline(kn * length_x / quantum_number_x, color="deepskyblue", linestyle=":", linewidth=1.4, alpha=0.85)

    if quantum_number_y > 1:
        has_2d_nodes = True
        for kn in range(1, quantum_number_y):
            ax.axhline(kn * length_y / quantum_number_y, color="deepskyblue", linestyle=":", linewidth=1.4, alpha=0.85)

    if has_2d_nodes:
        ax.plot([], [], color="deepskyblue", linestyle=":", linewidth=1.4, label="Nodal lines (rho = 0)")
        ax.legend(loc="upper right")

    ax.set_xlabel("Position x (m)")
    ax.set_ylabel("Position y (m)")
    ax.set_title(
        "2D Probability Density (nx = " + str(quantum_number_x) +
        ", ny = " + str(quantum_number_y) + ")"
    )
    plt.tight_layout()
    plt.show()


# ==============================================================================
# 3D WAVEFUNCTION ANALYSIS
# ==============================================================================

def wavefunction_3d():
    """Perform complete 3D potential well wavefunction analysis and transparent 3D plotting."""
    print("\n===================================")
    print("     3D WAVEFUNCTION ANALYSIS")
    print("===================================")
    print("Given 3D wavefunction:")
    print("psi(x, y, z) = A * sin(nx*pi*x/Lx) * sin(ny*pi*y/Ly) * sin(nz*pi*z/Lz)")
    print("where A = sqrt(8 / (Lx * Ly * Lz))")

    m = mass()
    length_x = length("the well length Lx")
    length_y = length("the well length Ly")
    length_z = length("the well length Lz")

    print("\nEnter quantum numbers (nx, ny, nz):")
    quantum_number_x = ask_int("Enter quantum number nx (1, 2, 3, ...): ")
    quantum_number_y = ask_int("Enter quantum number ny (1, 2, 3, ...): ")
    quantum_number_z = ask_int("Enter quantum number nz (1, 2, 3, ...): ")

    # Normalization constant: A = sqrt(8 / (Lx * Ly * Lz))
    A = math.sqrt(8.0 / (length_x * length_y * length_z))
    print("\nNormalization constant A =", format(A, ".6e"))

    # Energy calculation: E = (h^2 / (8*m)) * (nx^2/Lx^2 + ny^2/Ly^2 + nz^2/Lz^2)
    print("\n---------- PARTICLE ENERGY ----------")
    energy_val = energy(
        (quantum_number_x, quantum_number_y, quantum_number_z),
        m,
        (length_x, length_y, length_z)
    )
    show_energy(energy_val)

    # Probability density at a point (x, y, z)
    print("\n---------- PROBABILITY DENSITY AT A POINT ----------")
    while True:
        x = length("the position x (0 <= x <= Lx)", allow_zero=True)
        if x <= length_x or math.isclose(x, length_x, rel_tol=1e-12):
            break
        print("Position x must be between 0 and Lx.")

    while True:
        y = length("the position y (0 <= y <= Ly)", allow_zero=True)
        if y <= length_y or math.isclose(y, length_y, rel_tol=1e-12):
            break
        print("Position y must be between 0 and Ly.")

    while True:
        z = length("the position z (0 <= z <= Lz)", allow_zero=True)
        if z <= length_z or math.isclose(z, length_z, rel_tol=1e-12):
            break
        print("Position z must be between 0 and Lz.")

    psi = psi_value_3d(
        x, y, z,
        length_x, length_y, length_z,
        quantum_number_x, quantum_number_y, quantum_number_z,
        A
    )
    probability_density = psi**2

    print("Wavefunction psi(x, y, z) =", format(psi, ".6e"))
    print("Probability density =", format(probability_density, ".6e"), "per m^3")

    # Probability in a region [ax*Lx, bx*Lx] x [ay*Ly, by*Ly] x [az*Lz, bz*Lz]
    print("\n---------- PROBABILITY IN A REGION ----------")
    print("Enter limits along x as fractions of Lx (0 to 1).")
    while True:
        ax = ask_float("Enter lower limit ax/Lx (0 to 1): ")
        bx = ask_float("Enter upper limit bx/Lx (0 to 1): ")
        if 0 <= ax < bx <= 1:
            break
        print("Invalid limits. Require 0 <= ax/Lx < bx/Lx <= 1.")

    print("\nEnter limits along y as fractions of Ly (0 to 1).")
    while True:
        ay = ask_float("Enter lower limit ay/Ly (0 to 1): ")
        by = ask_float("Enter upper limit by/Ly (0 to 1): ")
        if 0 <= ay < by <= 1:
            break
        print("Invalid limits. Require 0 <= ay/Ly < by/Ly <= 1.")

    print("\nEnter limits along z as fractions of Lz (0 to 1).")
    while True:
        az = ask_float("Enter lower limit az/Lz (0 to 1): ")
        bz = ask_float("Enter upper limit bz/Lz (0 to 1): ")
        if 0 <= az < bz <= 1:
            break
        print("Invalid limits. Require 0 <= az/Lz < bz/Lz <= 1.")

    total_prob = probability_between_3d(
        quantum_number_x, quantum_number_y, quantum_number_z,
        ax, bx, ay, by, az, bz
    )
    print("Probability =", format(total_prob, ".6f"))
    print("Probability percentage =", format(total_prob * 100.0, ".4f"), "%")

    # 3D Probability density graph
    print("\n---------- 3D PROBABILITY DENSITY GRAPH ----------")

    # FIX: skip the graph when any n is too large to draw reliably
    if max(quantum_number_x, quantum_number_y, quantum_number_z) > PLOT_MAX_N_3D:
        print("Quantum numbers are too large to draw a reliable 3D graph (limit n <=", PLOT_MAX_N_3D, ").")
        print("Graph skipped. The numerical results above are still valid.")
        return

    print("Generating 3D graph (including zero-probability regions and nodal planes)...")

    # Balanced grid resolution for smooth visualization and responsive performance
    grid_points = min(max(22, 7 * max(quantum_number_x, quantum_number_y, quantum_number_z)), 35)
    rho_max = A**2

    x_points = []
    y_points = []
    z_points = []
    density_values = []

    for i in range(grid_points):
        x_val = length_x * (i + 0.5) / grid_points
        for j in range(grid_points):
            y_val = length_y * (j + 0.5) / grid_points
            for k in range(grid_points):
                z_val = length_z * (k + 0.5) / grid_points

                psi_pt = psi_value_3d(
                    x_val, y_val, z_val,
                    length_x, length_y, length_z,
                    quantum_number_x, quantum_number_y, quantum_number_z,
                    A
                )
                rho_pt = psi_pt**2

                x_points.append(x_val)
                y_points.append(y_val)
                z_points.append(z_val)
                density_values.append(rho_pt)

    fig = plt.figure(figsize=(9, 7))
    ax = fig.add_subplot(projection="3d")

    cmap = plt.get_cmap("plasma")   # FIX: works on all Matplotlib versions
    point_colors = []
    point_sizes = []

    # Map colors and transparency continuously across all points:
    # Zero/low density points are faint, small, and highly transparent (alpha ~ 0.03, s ~ 6),
    # while peak lobes are vivid and brighter (alpha ~ 0.65, s ~ 28).
    # This represents the 0-probability regions without forming an opaque, solid cuboid.
    for rho in density_values:
        norm_rho = min(1.0, max(0.0, rho / rho_max))
        rgba = list(cmap(norm_rho))
        rgba[3] = 0.03 + 0.62 * (norm_rho ** 0.7)
        point_colors.append(tuple(rgba))
        point_sizes.append(6.0 + 22.0 * norm_rho)

    sc = ax.scatter(
        x_points, y_points, z_points,
        color=point_colors,
        s=point_sizes,
        edgecolors="none"
    )

    # Wireframe bounding box represents the physical well boundaries where rho = 0 at the walls
    for x_edge in [0, length_x]:
        for y_edge in [0, length_y]:
            ax.plot([x_edge, x_edge], [y_edge, y_edge], [0, length_z], color="gray", linestyle="--", linewidth=0.8, alpha=0.5)
    for x_edge in [0, length_x]:
        for z_edge in [0, length_z]:
            ax.plot([x_edge, x_edge], [0, length_y], [z_edge, z_edge], color="gray", linestyle="--", linewidth=0.8, alpha=0.5)
    for y_edge in [0, length_y]:
        for z_edge in [0, length_z]:
            ax.plot([0, length_x], [y_edge, y_edge], [z_edge, z_edge], color="gray", linestyle="--", linewidth=0.8, alpha=0.5)

    # Delineate internal nodal planes where probability density is identically zero (rho = 0)
    has_nodal_planes = False
    if quantum_number_x > 1:
        has_nodal_planes = True
        for kn in range(1, quantum_number_x):
            x_node = kn * length_x / quantum_number_x
            for y_line in [0, 0.5 * length_y, length_y]:
                ax.plot([x_node, x_node], [y_line, y_line], [0, length_z], color="deepskyblue", linestyle=":", linewidth=1.0, alpha=0.6)
            for z_line in [0, 0.5 * length_z, length_z]:
                ax.plot([x_node, x_node], [0, length_y], [z_line, z_line], color="deepskyblue", linestyle=":", linewidth=1.0, alpha=0.6)

    if quantum_number_y > 1:
        has_nodal_planes = True
        for kn in range(1, quantum_number_y):
            y_node = kn * length_y / quantum_number_y
            for x_line in [0, 0.5 * length_x, length_x]:
                ax.plot([x_line, x_line], [y_node, y_node], [0, length_z], color="deepskyblue", linestyle=":", linewidth=1.0, alpha=0.6)
            for z_line in [0, 0.5 * length_z, length_z]:
                ax.plot([0, length_x], [y_node, y_node], [z_line, z_line], color="deepskyblue", linestyle=":", linewidth=1.0, alpha=0.6)

    if quantum_number_z > 1:
        has_nodal_planes = True
        for kn in range(1, quantum_number_z):
            z_node = kn * length_z / quantum_number_z
            for x_line in [0, 0.5 * length_x, length_x]:
                ax.plot([x_line, x_line], [0, length_y], [z_node, z_node], color="deepskyblue", linestyle=":", linewidth=1.0, alpha=0.6)
            for y_line in [0, 0.5 * length_y, length_y]:
                ax.plot([0, length_x], [y_line, y_line], [z_node, z_node], color="deepskyblue", linestyle=":", linewidth=1.0, alpha=0.6)

    # Legend clearly explaining physical zero and peak features
    ax.plot([], [], color="gray", linestyle="--", linewidth=1.0, label="Well boundary (rho = 0 at walls)")
    if has_nodal_planes:
        ax.plot([], [], color="deepskyblue", linestyle=":", linewidth=1.2, label="Nodal plane(s) (rho = 0)")
    ax.scatter([], [], color=cmap(0.0), s=15, alpha=0.6, label="Zero / minimum density")
    ax.scatter([], [], color=cmap(1.0), s=25, alpha=0.9, label="Peak probability density")
    ax.legend(loc="upper right", fontsize=8)

    # Colorbar starting at 0
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0.0, vmax=rho_max))
    sm.set_array([])
    cb = fig.colorbar(sm, ax=ax, pad=0.1)
    cb.set_label("Probability density |psi(x, y, z)|^2 (m^-3)")

    ax.set_xlabel("Position x (m)")
    ax.set_ylabel("Position y (m)")
    ax.set_zlabel("Position z (m)")
    ax.set_title(
        "3D Probability Density (nx = " + str(quantum_number_x) +
        ", ny = " + str(quantum_number_y) +
        ", nz = " + str(quantum_number_z) + ")"
    )

    plt.tight_layout()
    plt.show()


# ==============================================================================
# WAVEFUNCTION ANALYSIS ROUTER
# ==============================================================================

def wavefunction():
    """Wavefunction Analysis dimension selection menu."""
    print("\n===================================")
    print("       WAVEFUNCTION ANALYSIS")
    print("===================================")
    print("Select the dimension of the potential well:")
    print("1. One-dimensional (1D)")
    print("2. Two-dimensional (2D)")
    print("3. Three-dimensional (3D)")

    while True:
        choice = input("Enter your choice (1-3): ").strip()

        if choice == "1":
            wavefunction_1d()
            break
        if choice == "2":
            wavefunction_2d()
            break
        if choice == "3":
            wavefunction_3d()
            break

        print("Invalid dimension selection. Please enter 1, 2 or 3.")


# ==============================================================================
# POTENTIAL WELL TRANSITIONS (1D, 2D, and 3D)
# ==============================================================================

def transitions():
    """
    Calculate transition energies, absorption/emission classification,
    and associated photon properties in 1D, 2D, or 3D potential wells.
    """
    print("\n===================================")
    print("     POTENTIAL WELL TRANSITIONS")
    print("===================================")
    print("Select the dimension of the potential well:")
    print("1. One-dimensional (1D)")
    print("2. Two-dimensional (2D)")
    print("3. Three-dimensional (3D)")

    while True:
        choice = input("Enter your choice (1-3): ").strip()

        if choice == "1":
            d = 1
            break
        if choice == "2":
            d = 2
            break
        if choice == "3":
            d = 3
            break

        print("Invalid dimension selection. Please enter 1, 2 or 3.")

    m = mass()

    L = []
    initial = []
    final = []

    axes = ["x", "y", "z"]

    for i in range(d):
        L.append(length("the well length L" + axes[i]))

        print("\nEnter the INITIAL quantum number for", axes[i])
        initial.append(ask_int("Initial quantum number n" + axes[i] + ": "))

        print("\nEnter the FINAL quantum number for", axes[i])
        final.append(ask_int("Final quantum number n" + axes[i] + ": "))

    # Calculate energies
    Ei = energy(initial, m, L)
    Ef = energy(final, m, L)
    delta = energy_difference(initial, final, m, L)

    print("\n===================================")
    print("          TRANSITION RESULTS")
    print("===================================")
    print("Dimension:", str(d) + "D")
    print("Initial state:", tuple(initial))
    print("Final state:", tuple(final))

    print("\n---------- INITIAL ENERGY ----------")
    show_energy(Ei)

    print("\n---------- FINAL ENERGY ----------")
    show_energy(Ef)

    print("\n---------- SIGNED ENERGY DIFFERENCE ----------")
    show_energy(delta)

    # FIX: identical states are not "degenerate states", so check this first
    if initial == final:
        print("\nInitial and final states are identical. No transition occurs.")
        return

    if delta > 0:
        print("\nTransition: ABSORPTION")
        print("The particle absorbs energy.")
    elif delta < 0:
        print("\nTransition: EMISSION")
        print("The particle emits energy.")
    else:
        print("\nNo energy change occurs (the two different states are degenerate).")
        return

    print("\n---------- ENERGY EXCHANGED ----------")
    show_energy(abs(delta))

    # Photon properties
    frequency = abs(delta) / h
    wavelength = c / frequency

    print("\n---------- PHOTON INFORMATION ----------")
    print("Photon frequency =", format(frequency, ".6e"), "Hz")
    print("Photon wavelength =", format(wavelength, ".6e"), "m")
    print("Photon wavelength =", format(wavelength * 1e9, ".6f"), "nm")

    # FIX: reminder that selection rules are not applied
    print("\nNote: selection rules are not applied. For a charged particle in a")
    print("1D well, only transitions with an odd change in n are optically allowed.")


# ==============================================================================
# MAIN MENU
# ==============================================================================

def main():
    """Main application loop."""
    while True:
        print("\n======================================")
        print("      QUANTUM PARTICLE CALCULATOR")
        print("======================================")
        print("1. Wavefunction Analysis")
        print("2. Potential Well Transitions")
        print("3. Exit the Program")

        choice = input("\nEnter your choice (1-3): ").strip()

        if choice == "1":
            wavefunction()
        elif choice == "2":
            transitions()
        elif choice == "3":
            print("\nThank you for using the calculator.")
            print("Program terminated.")
            break
        else:
            print("Invalid choice. Please enter 1, 2 or 3.")


if __name__ == "__main__":
    # FIX: exit cleanly on Ctrl+C / Ctrl+D instead of a traceback
    try:
        main()
    except (EOFError, KeyboardInterrupt):
        print("\nProgram terminated.")
