import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
import io

# --- Page Setup ---
st.set_page_config(page_title="Cosmic Explorer", layout="wide")

# -- GENUINELY have no idea why this isn't working, tabs aren't being updated. adjust spacing to liking and tix tabs 1 and 2. 3,4 and 5 haven't been added yet

tab1, tab2, tab3, tab4, tab5, tab6 = st.tabs([
    "Home",
    "Cosmology Calculator",
    "Galaxy Size",
    "Graphs",
    "Export",
    "More Tools"
])

# --- Add Banner Image (top-of-page full-width Hubble image) ---

# ----

# --- Optional Logo (centered telescope) ---

# ----


# --- Title and Subtitle ---
with tab1:
    st.markdown("""
    # 🌌 **Cosmic Explorer Observatory Toolbox**
    ##### 🧠 Real-Time Cosmology Calculator · Galaxy Size Estimator · Distance Grapher
    ---
    """)
# --------------------------- ADD INFORMATIONAL PAGE AS PART OF HOME (TAB 1) --------------------------
# --- REMOVE WHEN FIXED: Banners and Logo don't work; remove or fix

# --- Sidebar Sliders ---
st.sidebar.header("🔧 Cosmology Inputs")
H0 = st.sidebar.slider("Hubble Constant (H₀)", 50, 90, 70,
                       help="Expansion rate of the universe today (in km/s/Mpc)")
Om0 = st.sidebar.slider("Matter Density (Ωₘ)", 0.0, 1.0, 0.3,
                        help="Fraction of total density due to matter")
Ol0 = st.sidebar.slider("Dark Energy Density (ΩΛ)", 0.0, 1.0, 0.7,
                        help="Fraction of total density due to dark energy")
z = st.sidebar.slider("Redshift (z)", 0.0, 10.0, 2.0,
                      help="How far back in time you’re looking (z=0=today)")


cosmo = FlatLambdaCDM(H0=H0 * u.km / u.s / u.Mpc, Om0=Om0)

# --- Cosmology Outputs ---
with tab2: 
    age = cosmo.age(0).value
    age_z = cosmo.age(z).value
    lookback = cosmo.lookback_time(z).value
    comoving = cosmo.comoving_distance(z).value
    angular = cosmo.angular_diameter_distance(z).value
    luminosity = cosmo.luminosity_distance(z).value
   
    st.header("🧮 Cosmology Outputs")
    st.markdown(f"- **Age of Universe now:** {age:.2f} Gyr")
    st.markdown(f"- **Age at z={z:.2f}:** {age_z:.2f} Gyr")
    st.markdown(f"- **Lookback Time:** {lookback:.2f} Gyr")
    st.markdown(f"- **Comoving Distance:** {comoving:.2f} Mpc")
    st.markdown(f"- **Angular Diameter Distance:** {angular:.2f} Mpc")
    st.markdown(f"- **Luminosity Distance:** {luminosity:.2f} Mpc")


# --- 📏 Galaxy Size Calculator ---
with tab3:
    st.markdown("---")
    st.header("📏 Galaxy Size Calculator (Angular → Physical Size)")

    # Input: Angular size in arcseconds
    theta_arcsec = st.number_input("Angular size (arcseconds)", min_value=0.0, value=30.0, step=1.0)

    # Input: Unit toggle
    unit_choice = st.radio("Display physical size in:", ("kpc", "light-years", "both"))

    # Calculation
    theta_rad = theta_arcsec * (np.pi / (180 * 3600))  # arcsec to radians
    size_kpc = angular * 1000 * theta_rad  # angular is in Mpc → *1000 to kpc
    size_ly = size_kpc * 3261.56  # 1 kpc ≈ 3261.56 light-years

    # Output
    st.subheader("🪐 Estimated Physical Size")

    if unit_choice == "kpc":
        st.markdown(f"- **Size**: {size_kpc:.2f} kiloparsecs")
    elif unit_choice == "light-years":
        st.markdown(f"- **Size**: {size_ly:.2f} light-years")
    else:
        st.markdown(f"- **Size**: {size_kpc:.2f} kpc ≈ {size_ly:.2f} light-years")

    # --- Save Galaxy Size Result ---
    import io

    # Create a simple text or CSV string
    if unit_choice == "kpc":
        size_result = f"{size_kpc:.2f} kpc"
    elif unit_choice == "light-years":
        size_result = f"{size_ly:.2f} light-years"
    else:
        size_result = f"{size_kpc:.2f} kpc ≈ {size_ly:.2f} light-years"

    save_text = f"""Galaxy Size Calculation
    -----------------------
    Redshift (z): {z:.2f}
    Angular Size: {theta_arcsec:.2f} arcseconds
    Physical Size: {size_result}
    """

    # Convert to bytes for download
    buffer = io.BytesIO()
    buffer.write(save_text.encode())
    buffer.seek(0)

    # Download button
    st.download_button(
        label="💾 Save Galaxy Size Result",
        data=buffer,
        file_name="galaxy_size.txt",
        mime="text/plain"
    )

# --- Lookback Time Graph ---
# --- Lookback Time Plot ---
with tab4:
    st.markdown("---")
    st.header("⏳ Lookback Time vs Redshift")

    z_range = np.linspace(0.01, 10, 500)
    lookback_range = cosmo.lookback_time(z_range).value

    fig1, ax1 = plt.subplots()
    ax1.plot(z_range, lookback_range)
    ax1.set_xlabel("Redshift (z)")
    ax1.set_ylabel("Lookback Time (Gyr)")
    ax1.grid(True)
    ax1.set_title("Lookback Time vs Redshift")
    st.pyplot(fig1)

    # Save PNG to buffer
    lookback_buf = io.BytesIO()
    fig1.savefig(lookback_buf, format="png")
    lookback_buf.seek(0)

    # Download Button
    st.download_button(
        label="💾 Save Lookback Time Graph as PNG",
        data=lookback_buf,
        file_name="lookback_time.png",
        mime="image/png"
    )


    # --- Distance Graph ---
    # --- Distance Plot ---
    st.markdown("---")
    st.header("📏 Distances vs Redshift")
    st.markdown("""
    Comparing the distance and redshift shows us the **how much light from distant galaxies has been stretched (z) due to cosmic expansion in relation to distance.**
    - Comoving Distance is the fixed present day distance to an object (ignores expansion after light left).
    - Angular Diameter Distance matches how an object's size appears on the sky.
    - Luminosity Distance is an inferred distance from how dim an object appears.
    """)
    com = cosmo.comoving_distance(z_range).value / 1000
    ang = cosmo.angular_diameter_distance(z_range).value / 1000
    lum = cosmo.luminosity_distance(z_range).value / 1000

    fig2, ax2 = plt.subplots()
    ax2.plot(z_range, com, label="Comoving (Glyr)")
    ax2.plot(z_range, ang, label="Angular Diameter (Glyr)")
    ax2.plot(z_range, lum, label="Luminosity (Glyr)")
    ax2.set_xlabel("Redshift (z)")
    ax2.set_ylabel("Distance (Glyr)")
    ax2.legend()
    ax2.grid(True)
    ax2.set_title("Distances vs Redshift")
    st.pyplot(fig2)

    # Save PNG to buffer
    distance_buf = io.BytesIO()
    fig2.savefig(distance_buf, format="png")
    distance_buf.seek(0)

    # Download Button
    st.download_button(
        label="💾 Save Distance Graph as PNG",
        data=distance_buf,
        file_name="distances_vs_redshift.png",
        mime="image/png"
    )

     # 📈 Expansion History Plot
    st.subheader("📈 Expansion History")
    st.markdown("""
    This graph shows how the **scale factor** (size of the universe) has changed over time.
    - The universe expands faster as it gets older
    - Early on, expansion was slower due to gravity's influence
    """)
    z_exp = np.linspace(0.01, 10, 500)
    t_exp = cosmo.age(z_exp).value
    a_exp = 1 / (1 + z_exp)

    fig_exp, ax_exp = plt.subplots()
    sc = ax_exp.scatter(t_exp, a_exp, c=z_exp, cmap='plasma', s=10)
    ax_exp.set_xlabel("Time Since Big Bang (Gyr)")
    ax_exp.set_ylabel("Scale Factor (a)")
    ax_exp.set_title("Expansion History of the Universe")
    ax_exp.grid(True)
    cbar = plt.colorbar(sc, ax=ax_exp, label="Redshift (z)")
    st.pyplot(fig_exp)

    # 💾 Save Expansion History Graph
    expansion_buf = io.BytesIO()
    fig_exp.savefig(expansion_buf, format="png")
    expansion_buf.seek(0)

    st.download_button(
        label="💾 Download Expansion History Graph (PNG)",
        data=expansion_buf,
        file_name="expansion_history.png",
        mime="image/png"
    )

# --- 📄 Summary Report Exporter ---
with tab5:
    st.markdown("---")
    st.header("📄 Export Summary Report")

    # Recompute galaxy size in both units
    theta_rad = theta_arcsec * (np.pi / (180 * 3600))
    size_kpc = angular * 1000 * theta_rad
    size_ly = size_kpc * 3261.56

    # Generate report text
    report = f"""
    Cosmic Explorer Summary Report
    ------------------------------
    Cosmology Inputs:
    - Hubble Constant (H₀): {H0} km/s/Mpc
    - Matter Density (Ωₘ): {Om0}
    - Dark Energy Density (ΩΛ): {Ol0}
    - Redshift (z): {z}

    Cosmology Results:
    - Current Age of Universe: {age:.2f} Gyr
    - Age at z = {z:.2f}: {age_z:.2f} Gyr
    - Lookback Time: {lookback:.2f} Gyr
    - Comoving Distance: {comoving:.2f} Mpc
    - Angular Diameter Distance: {angular:.2f} Mpc
    - Luminosity Distance: {luminosity:.2f} Mpc

    Galaxy Size Calculation:
    - Angular Size: {theta_arcsec:.2f} arcseconds
    - Physical Size: {size_kpc:.2f} kpc ≈ {size_ly:.2f} light-years
    """

    # Convert to downloadable buffer
    report_buf = io.BytesIO()
    report_buf.write(report.encode())
    report_buf.seek(0)

    # Download button
    st.download_button(
        label="📄 Download Summary Report",
        data=report_buf,
        file_name="cosmic_summary.txt",
        mime="text/plain"
    )

with tab6:
    st.header("🧪 More Tools")
    st.subheader("📅 Cosmic Time at Redshift")

    # Age of universe at z is already computed!
    st.write(f"At redshift z = {z:.2f}, the universe was approximately **{age_z:.2f} Gyr** old.")

    # 🔭 Independent redshift slider
    st.markdown("Adjust the redshift independently of the main app:")
    z_tools = st.slider("Select redshift (z)", 0.01, 10.0, 2.0, step=0.01)

    # 📉 Scale Factor
    st.subheader("Scale Factor (a)")
    st.markdown("""
    The **scale factor** `a = 1 / (1 + z)` tells us how much the universe has expanded.
    - At `a = 1`, the universe is at its current size.
    - At `a = 0.5`, the universe was half as large as it is now.
    """)
    a_tools = 1 / (1 + z_tools)
    st.write(f"At redshift **z = {z_tools:.2f}**, the scale factor was **a = {a_tools:.4f}**")

    # ⌛ Time Since Big Bang
    st.subheader("Time Since the Big Bang")
    st.markdown("""
    This is how much time had passed since the Big Bang at a given redshift.
    It decreases as you go further back in time (higher `z`).
    """)
    age_at_z_tools = cosmo.age(z_tools).value
    st.write(f"At redshift **z = {z_tools:.2f}**, the universe was approximately **{age_at_z_tools:.2f} Gyr** old.")

    # 🌌 Universe Timeline Summary
    st.subheader("Universe Timeline Summary")
    st.markdown("""
    | Time After Big Bang | Event |
    |---------------------|------------------------------|
    | ~0 s                | Big Bang Singularity         |
    | ~10⁻³⁵ s            | Inflation begins             |
    | ~10⁻⁶ s             | Quarks form protons/neutrons |
    | ~3 minutes          | Nucleosynthesis (H, He)      |
    | ~380,000 years      | Recombination / CMB released |
    | ~100 million years  | First stars form             |
    | ~1 billion years    | First galaxies               |
    | ~13.8 billion years | Today                        |
    """)

    # 🚀 Light Travel Distance
    st.subheader("Light Travel Distance")
    st.markdown("""
    This is the distance light has traveled from an object at redshift `z` to reach us today.
    It's often slightly **less** than the comoving distance.
    """)

    light_travel_time = cosmo.lookback_time(z_tools).value  # in Gyr
    light_travel_distance_mpc = light_travel_time * 0.306601  # 1 Gyr ~ 0.306601 Mpc (speed of light)

    st.write(f"At redshift **z = {z_tools:.2f}**, light has traveled approximately:")
    st.write(f"- **{light_travel_time:.2f} billion years (Gyr)**")
    st.write(f"- **{light_travel_distance_mpc:.2f} megaparsecs (Mpc)**")

with tab1:
    st.markdown("##### Made with ❤️ by Bri · Powered by Streamlit + Astropy")

