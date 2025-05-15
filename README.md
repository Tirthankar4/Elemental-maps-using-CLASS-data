# Lunar Elemental Mapper: Chandrayaan-2 CLASS Data Analysis  
*High-resolution mapping of lunar surface composition using X-ray fluorescence (XRF) data*

---

## 🚀 Project Overview  
This repository automates the analysis of **Chandrayaan-2 CLASS L1 data** to generate high-resolution elemental ratio maps (e.g., Mg/Si, Al/Si) for lunar geology studies. The pipeline identifies solar flare-associated XRF signals, processes spectral data, and creates sub-pixel resolution maps (4.125 km²) using a custom Drizzle-inspired algorithm.

**Key Scientific Goal**  
Map elemental heterogeneity to study lunar crust evolution and volcanic history through ratios of Mg, Al, Ca, Fe, Ti, Na, and O relative to Si.

---

## 🔧 Technical Components  
### Core Scripts  
| Script | Function |  
|--------|----------|  
| `fetch-flares.py` | Extracts CLASS L1 data matching flare windows (+15min buffer) |  
| `addClassSpectra.py` | Aggregates 8s spectra into time-integrated FITS files |  
| `gen_bg.py` | Generates OGIP-compliant background files from nighttime data |  
| `catalogue.py` | Automated peak detection via Gaussian Mixture Models (R² > 0.85) |  
| `colorMap.py` | Generates publication-ready elemental ratio maps in equirectangular projection |  

### Data Pipeline  
1. **Flare Sorting**: Uses Suryadrishti/Hinode XRT catalogs to isolate flare periods  
2. **SPICE Processing**: Identifies lunar day/night phases for background subtraction  
3. **Spectral Correction**: Applies ARF/RMF matrices via PyXSPEC  
4. **Sub-Pixel Mapping**: Implements custom Drizzle algorithm (4x resolution boost)  

---

## 📊 Key Outputs  
- **Elemental Ratio Maps**  
  ![Mg/Si Coverage Example][https://drive.google.com/uc?id=1tCE5y6DMkdBzSJCiCnDHJzPflYCuIkcm]
  *Sub-pixel resolution (4.125 km) maps for 7 elements across X/M/C-class flares*

- **Interactive Visualization**  
  - 2D/3D web-based viewer with layer toggling (Al/Si, Mg/Si, etc.)  
  - Coordinate-aware exploration of lunar maria/highlands  

---


**Dependencies**  
- HEASoft + PyXSPEC for spectral analysis  
- Astropy for FITS handling  
- GeoPandas + Cartopy for spatial analysis  

---

## 🧩 Challenges Addressed  
- **Data Volume**: Processed 5TB+ CLASS L1 data (2019-2024)  
- **ARF/RMF Corrections**: Fixed OGIP header compliance issues  
- **QGIS Limitations**: Switched to matplotlib/Blender for CRS-free mapping  

---


---

## 🔮 Future Work  
- Machine learning-based noise reduction  
- API integration with ISSDC Pradan portal  
- Augmented Reality visualization module  

*Developed for Inter-IIT Tech Meet 13.0 (Team 62) in collaboration with ISRO*

