import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
import io

# --- Page Setup ---
st.set_page_config(page_title="Cosmic Explorer", layout="wide")

# --- Add Banner Image (top-of-page full-width Hubble image) ---


# --- Optional Logo (centered telescope) ---




# --- Title and Subtitle ---
st.markdown("""
# 🌌 **Cosmic Explorer Observatory Toolbox**
##### 🧠 Real-Time Cosmology Calculator · Galaxy Size Estimator · Distance Grapher
---
""")

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
st.subheader("🧮 Cosmology Outputs")
age = cosmo.age(0).value
age_z = cosmo.age(z).value
lookback = cosmo.lookback_time(z).value
comoving = cosmo.comoving_distance(z).value
angular = cosmo.angular_diameter_distance(z).value
luminosity = cosmo.luminosity_distance(z).value

st.markdown(f"- **Age of Universe now:** {age:.2f} Gyr")
st.markdown(f"- **Age at z={z:.2f}:** {age_z:.2f} Gyr")
st.markdown(f"- **Lookback Time:** {lookback:.2f} Gyr")
st.markdown(f"- **Comoving Distance:** {comoving:.2f} Mpc")
st.markdown(f"- **Angular Diameter Distance:** {angular:.2f} Mpc")
st.markdown(f"- **Luminosity Distance:** {luminosity:.2f} Mpc")


# --- 📏 Galaxy Size Calculator ---
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

# --- 📄 Summary Report Exporter ---
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

st.markdown("##### Made with ❤️ by Bri · Powered by Streamlit + Astropy")

