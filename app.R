# app.R
# Nakad et al. (2023): Stomatal optimality, VPD response, and AFM (exploratory tool)
#
# This version:
#  - Keeps the simplified Nakad framework for VPD-response curves + AFM diagnostics.
#  - Reports "true" marginal WUE as in Nakad/Katul:
#      mWUE ≡ (∂f_c/∂f_e)|D = (∂f_c/∂g)/(∂f_e/∂g) = λ
#    (shown as text + consistency check; not plotted vs VPD because it is constant here).
#  - Adds a new tab: "A–Ci (Farquhar)" showing Ac (Rubisco limitation), Aj (RuBP limitation),
#    and A = min(Ac,Aj) - Rd, with an operating Ci marker from the Nakad-optimality Ci/Ca at a chosen VPD.

library(shiny)
library(ggplot2)

# -------------------------
# Helpers
# -------------------------
clamp <- function(x, lo, hi) pmin(pmax(x, lo), hi)

vline_if_in_range <- function(x, xmin, xmax, linetype = "dashed") {
  if (!is.na(x) && is.finite(x) && x >= xmin && x <= xmax) {
    geom_vline(xintercept = x, linetype = linetype)
  } else {
    NULL
  }
}

# numeric derivative dy/dx (finite differences), returns same length as x (first NA)
dydx_fd <- function(x, y) {
  dx <- diff(x)
  dy <- diff(y)
  c(NA_real_, dy / dx)
}

# -------------------------
# Model functions (simplified Nakad)
# -------------------------

# Eq. 12: optimal stomatal conductance to CO2 (g_CO2)
g_opt <- function(D_molmol, ca_molmol, a, lambda, alpha1_mol, alpha2_molmol, s) {
  D_molmol <- pmax(D_molmol, 1e-12)
  lambda <- pmax(lambda, 1e-12)
  
  g <- (alpha1_mol / (alpha2_molmol + s * ca_molmol)) *
    (-1 + sqrt(ca_molmol / (a * lambda * D_molmol)))
  
  pmax(g, 0)
}

# Eq. 10 (linearized demand): photosynthesis fc(g) in mol CO2 m^-2 s^-1
fc_from_g <- function(g_mol, ca_molmol, alpha1_mol, alpha2_molmol, s) {
  (g_mol * alpha1_mol * ca_molmol) / (alpha1_mol + g_mol * (alpha2_molmol + s * ca_molmol))
}

# Eq. 9: transpiration proxy fe ≈ a g D in mol H2O m^-2 s^-1
fe_from_g <- function(g_mol, D_molmol, a) {
  a * g_mol * D_molmol
}

# Eq. 20: optimal ci/ca (dimensionless) -- valid in open-stomata regime
ci_ca_opt <- function(D_molmol, ca_molmol, a, lambda) {
  D_molmol <- pmax(D_molmol, 1e-12)
  lambda <- pmax(lambda, 1e-12)
  1 - sqrt(a * lambda / ca_molmol) * sqrt(D_molmol)
}

# Invert Eq. 20 exactly to infer lambda from (ci/ca at Dref, ca):
# => lambda = (ca/(a*D))*(1 - ci/ca)^2
lambda_from_ci_ca_eq20 <- function(ci_ca, D_molmol, ca_molmol, a) {
  D_molmol <- pmax(D_molmol, 1e-12)
  ci_ca <- pmin(pmax(ci_ca, 0.05), 0.999)
  (ca_molmol / (a * D_molmol)) * (1 - ci_ca)^2
}

# Eq. 21: WUE (mol CO2 / mol H2O)
wue_closed_eq21 <- function(D_molmol, ca_molmol, a, lambda) {
  D_molmol <- pmax(D_molmol, 1e-12)
  lambda <- pmax(lambda, 1e-12)
  (ca_molmol / sqrt(a * D_molmol)) * sqrt(lambda / ca_molmol)
}

# Simplified AFM onset expression (with D in mol/mol):
# Dafm = (4/9) (ca/a) lambda^-1   [in D units]
D_afm_molmol <- function(ca_molmol, a, lambda) {
  lambda <- pmax(lambda, 1e-12)
  (4/9) * (ca_molmol / a) * (1 / lambda)
}

# Eq. 16-like sign term
Ca_afm_coeff <- function(alpha2_prime_molmol, ca_molmol, cp_molmol) {
  - (alpha2_prime_molmol - 2 * ca_molmol + 2 * cp_molmol) / (alpha2_prime_molmol^2)
}

# -------------------------
# TRUE marginal WUE (Nakad/Katul): ∂fc/∂fe at fixed D
# -------------------------
# For fc = (g α1 ca)/(α1 + g B), B = α2 + s ca:
# ∂fc/∂g = (α1^2 ca)/(α1 + gB)^2
dfc_dg <- function(g_mol, ca_molmol, alpha1_mol, alpha2_molmol, s) {
  B <- alpha2_molmol + s * ca_molmol
  denom <- alpha1_mol + g_mol * B
  (alpha1_mol^2 * ca_molmol) / (denom^2)
}

# fe = a g D => ∂fe/∂g = a D
dfe_dg <- function(D_molmol, a) {
  a * D_molmol
}

# -------------------------
# Optional: Farquhar A–Ci (two limitations)
# -------------------------
# Inputs/outputs in ppm & µmol m^-2 s^-1 for plotting.
farquhar_aci <- function(Ci_ppm,
                         Vcmax, J, Rd,
                         GammaStar_ppm,
                         Kc_ppm, Ko_ppm, O2_ppm) {
  
  # Rubisco-limited (Ac)
  Ac <- Vcmax * (Ci_ppm - GammaStar_ppm) / (Ci_ppm + Kc_ppm * (1 + O2_ppm / Ko_ppm))
  
  # RuBP-regeneration limited (Aj) common form
  Aj <- J * (Ci_ppm - GammaStar_ppm) / (4 * Ci_ppm + 8 * GammaStar_ppm)
  
  A <- pmin(Ac, Aj) - Rd
  
  data.frame(Ci_ppm = Ci_ppm, Ac = Ac, Aj = Aj, A = A)
}

# -------------------------
# Nakad Table 1 Vcmax(T)
# -------------------------
vcmax_T_table1 <- function(T_C, Vcmax25, m1, m2) {
  Vcmax25 * exp(m1 * (T_C - 25)) *
    (1 + exp(m2 * (T_C - 41)))^(-1)
}

# -------------------------
# λ(REW) (drought)
# -------------------------
lambda_from_REW <- function(REW, lambda_wet, form = c("power","linear"), k = 1.5, beta = 2) {
  form <- match.arg(form)
  REW <- pmax(pmin(REW, 1), 1e-6)
  
  if (form == "power") {
    # λ = λwet * REW^-k
    lambda_wet * (REW^(-k))
  } else {
    # λ = λwet * (1 + β (1-REW))
    lambda_wet * (1 + beta * (1 - REW))
  }
}

# -------------------------
# UI
# -------------------------
ui <- fluidPage(
  
  # Sticky sidebar CSS
  tags$head(
    tags$style(HTML("
      #sidebar {
        position: sticky;
        top: 10px;
        max-height: calc(100vh - 20px);
        overflow-y: auto;
      }
    "))
  ),
  
  titlePanel("Nakad et al. (2023): Stomatal optimality, VPD response, and AFM (exploratory tool)"),
  
  sidebarLayout(
    sidebarPanel(
      id = "sidebar",
      h4("Inputs"),
      
      sliderInput("ca_ppm", "Atmospheric CO₂, ca (ppm)", min = 250, max = 1000, value = 420, step = 10),
      
      checkboxInput("infer_lambda", "Infer λ from ci/ca (diagnostic mode; Eq. 20 inversion)", value = TRUE),
      
      conditionalPanel(
        condition = "input.infer_lambda == false",
        sliderInput("lambda_in", "λ (input; Lagrange multiplier; simplified units)",
                    min = 1e-5, max = 0.5, value = 0.03, step = 0.001)
      ),
      
      conditionalPanel(
        condition = "input.infer_lambda == true",
        sliderInput("ci_ca_ref", HTML("Reference c<sub>i</sub>/c<sub>a</sub> at VPD<sub>ref</sub>"),
                    min = 0.2, max = 0.95, value = 0.7, step = 0.01),
        sliderInput("VPD_ref", HTML("VPD<sub>ref</sub> (kPa) used to infer &lambda;"),
                    min = 0.1, max = 5, value = 1.0, step = 0.1)
      ),
      
      conditionalPanel(
        condition = "input.use_T_for_alpha1 == false",
        sliderInput(
          "alpha1_umol",
          label = HTML("&alpha;<sub>1</sub> (&mu;mol CO<sub>2</sub> m<sup>-2</sup> s<sup>-1</sup>; Vcmax-like)"),
          min = 10, max = 200, value = 80, step = 1
        )
      ),
      conditionalPanel(
        condition = "input.use_T_for_alpha1 == true",
        helpText("α1 is computed from Vcmax(T) (Table 1). To use the α1 slider, uncheck the box below.")
      ),
      
      sliderInput(
        "alpha2_ppm",
        label = HTML("&alpha;<sub>2</sub> (ppm; RuBisCO term like K<sub>c</sub>(1+O/K<sub>o</sub>) under RuBisCO limitation)"),
        min = 50, max = 1200, value = 550, step = 10
      ),
      
      sliderInput("s", "s (ci/ca linearization parameter)", min = 0.3, max = 0.9, value = 0.7, step = 0.01),
      
      sliderInput("a",
                  HTML("a (diffusivity ratio D<sub>H2O</sub>/D<sub>CO2</sub>; used in f<sub>e</sub> = a·g<sub>CO2</sub>·D)"),
                  min = 1.2, max = 2.0, value = 1.6, step = 0.01),
      
      sliderInput("patm_kPa", "Atmospheric pressure, Patm (kPa)", min = 80, max = 110, value = 101.3, step = 0.1),
      
      hr(),
      h4("VPD range"),
      sliderInput("VPDmin", "VPD min (kPa)", min = 0, max = 5, value = 0.1, step = 0.05),
      sliderInput("VPDmax", "VPD max (kPa)", min = 0.2, max = 10, value = 5, step = 0.1),
      numericInput("npts", "Number of points", value = 500, min = 100, max = 2000, step = 50),
      
      hr(),
      h4("Treatment presets"),
      selectInput(
        "preset",
        "Treatment preset",
        choices = c(
          "Custom",
          "Ambient (A)",
          "Heat (H)",
          "Drought (D)",
          "Heat + Drought (H+D)"
        ),
        selected = "Ambient"
      ),
      
      hr(),
      h4("Drought (soil water)"),
      sliderInput("REW", "Relative extractable water (REW)", min = 0.05, max = 1, value = 1, step = 0.01),
      selectInput("drought_fun", "λ(REW) form", choices = c("Power law" = "power", "Linear" = "linear"), selected = "power"),
      numericInput("lambda_wet", "λ at REW = 1 (wet reference)", value = 0.008, min = 1e-5, max = 1, step = 0.001),
      conditionalPanel(
        condition = "input.drought_fun == 'power'",
        numericInput("k_rew", "Drought sensitivity k (larger = stronger response)", value = 1.5, min = 0, max = 10, step = 0.1)
      ),
      conditionalPanel(
        condition = "input.drought_fun == 'linear'",
        numericInput("beta_rew", "Drought sensitivity β (λ = λwet*(1+β(1-REW)))", value = 2, min = 0, max = 20, step = 0.2)
      ),
      
      hr(),
      h4("Heat (temperature)"),
      sliderInput("Tair_C", "Air/leaf temperature T (°C)", min = 0, max = 50, value = 25, step = 0.5),
      checkboxInput("use_T_for_alpha1", "Use Table-1 Vcmax(T) for α1 (recommended)", value = TRUE),
      numericInput("Vcmax25", "Vcmax,25 (µmol m⁻² s⁻¹)", value = 80, min = 1, max = 300, step = 1),
      numericInput("m1", "m1 (°C⁻¹)", value = 0.08, min = 0, max = 0.3, step = 0.001),
      numericInput("m2", "m2 (°C⁻¹)", value = 0.20, min = 0, max = 1, step = 0.01),
      
      hr(),
      h4("Axis limits"),
      checkboxInput("fix_axes", "Fix y-axis limits (for comparisons)", value = TRUE),
      
      conditionalPanel(
        condition = "input.fix_axes == true",
        tags$div(
          tags$strong("g (CO₂) limits"),
          numericInput("gco2_min", "min (mol CO2 m^-2 s^-1)", value = 0, min = 0),
          numericInput("gco2_max", "max (mol CO2 m^-2 s^-1)", value = 0.5, min = 0),
          
          tags$hr(),
          
          tags$strong("g (H₂O) limits"),
          numericInput("gh2o_min", "min (mol H2O m^-2 s^-1)", value = 0, min = 0),
          numericInput("gh2o_max", "max (mol H2O m^-2 s^-1)", value = 0.8, min = 0),
          
          tags$hr(),
          
          tags$strong("f_c limits"),
          numericInput("fc_min", "min (µmol CO2 m^-2 s^-1)", value = 0, min = 0),
          numericInput("fc_max", "max (µmol CO2 m^-2 s^-1)", value = 50, min = 0),
          
          tags$hr(),
          
          tags$strong("f_e limits"),
          numericInput("fe_min", "min (mmol H2O m^-2 s^-1)", value = 0, min = 0),
          numericInput("fe_max", "max (mmol H2O m^-2 s^-1)", value = 5, min = 0),
          
          tags$hr(),
          
          tags$strong("c_i/c_a limits"),
          numericInput("cica_min", "min (-)", value = 0, min = 0, max = 1),
          numericInput("cica_max", "max (-)", value = 1, min = 0, max = 1),
          
          tags$hr(),
          
          tags$strong("WUE limits"),
          numericInput("wue_min", "min (mol/mol)", value = 0, min = 0),
          numericInput("wue_max", "max (mol/mol)", value = 0.02, min = 0)
        )
      ),
      
      hr(),
      h4("A–Ci (Farquhar) settings"),
      sliderInput("VPD_for_aci", "VPD for operating Ci marker (kPa)", min = 0.1, max = 5, value = 1.0, step = 0.1),
      sliderInput("Vcmax", "Vcmax (µmol m^-2 s^-1)", min = 10, max = 300, value = 80, step = 1),
      sliderInput("J", "J (or Jmax proxy) (µmol m^-2 s^-1)", min = 10, max = 400, value = 140, step = 5),
      sliderInput("Rd", "Rd (µmol m^-2 s^-1)", min = 0, max = 5, value = 1.0, step = 0.1),
      sliderInput("GammaStar", HTML("&Gamma;* (CO<sub>2</sub> compensation point, ppm)"), min = 10, max = 80, value = 42, step = 1),
      sliderInput("Kc", "Kc (ppm)", min = 50, max = 1500, value = 404, step = 10),
      sliderInput("Ko", "Ko (ppm)", min = 10000, max = 500000, value = 278000, step = 5000),
      sliderInput("O2", "O2 (ppm)", min = 150000, max = 230000, value = 210000, step = 1000),
      
      hr(),
      h4("Notes"),
      helpText(
        "We convert VPD (kPa) to D = VPD/Patm (dimensionless mol/mol) to keep Eq. 9 unit-consistent.",
        "Transpiration proxy: fe ≈ a g_CO2 D (where a≈1.6 accounts for diffusivity ratio).",
        "AFM: dashed line is theory Dafm; dotted line is numeric peak of fe within chosen VPD range.",
        "cp is used only for the Eq.16-like AFM sign diagnostic in this simplified implementation.",
        "Eq.20-based ci/ca is shown in the open-stomata regime; when g=0 we set ci/ca → 1 (physical limit).",
        "Marginal WUE in Nakad/Katul is defined as ∂fc/∂fe at fixed D, and equals λ."
      ),
      width = 4
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Curves",
                 fluidRow(
                   column(6, plotOutput("p_gco2", height = 300)),
                   column(6, plotOutput("p_gh2o", height = 300))
                 ),
                 fluidRow(
                   column(6, plotOutput("p_fc", height = 300)),
                   column(6, plotOutput("p_fe", height = 300))
                 ),
                 fluidRow(
                   column(6, plotOutput("p_ci", height = 300)),
                   column(6, plotOutput("p_wue", height = 300))
                 )
        ),
        tabPanel("AFM diagnostics",
                 plotOutput("p_dfe"),
                 verbatimTextOutput("txt_afm"),
                 hr(),
                 verbatimTextOutput("txt_mwue")
        ),
        tabPanel("A–Ci (Farquhar)",
                 plotOutput("p_aci", height = 450),
                 helpText("Ac: Rubisco limitation; Aj: RuBP regeneration limitation; A = min(Ac,Aj) - Rd. Dashed vertical line marks Ci predicted by Nakad-optimality at the chosen VPD.")
        ),
        tabPanel("WUE geometry",
                 plotOutput("p_wue_geom", height = 500),
                 helpText("Solid curve: fc vs fe as VPD varies (solid: pre-AFM branch, dotted: post-AFM branch). 
                  Dashed line: average slope (WUE = fc/fe at selected point). 
                  Solid straight line: tangent slope (mWUE = λ).")
        )
      ),
      width = 8
    )
  )
)

# -------------------------
# Server
# -------------------------
server <- function(input, output, session) {
  
  # y-limits (optional fixed axes)
  ylim_gco2 <- reactive({
    if (!isTRUE(input$fix_axes)) return(NULL)
    c(min(input$gco2_min, input$gco2_max), max(input$gco2_min, input$gco2_max))
  })
  ylim_gh2o <- reactive({
    if (!isTRUE(input$fix_axes)) return(NULL)
    c(min(input$gh2o_min, input$gh2o_max), max(input$gh2o_min, input$gh2o_max))
  })
  ylim_fc <- reactive({
    if (!isTRUE(input$fix_axes)) return(NULL)
    c(min(input$fc_min, input$fc_max), max(input$fc_min, input$fc_max))
  })
  ylim_fe <- reactive({
    if (!isTRUE(input$fix_axes)) return(NULL)
    c(min(input$fe_min, input$fe_max), max(input$fe_min, input$fe_max))
  })
  ylim_cica <- reactive({
    if (!isTRUE(input$fix_axes)) return(NULL)
    c(clamp(min(input$cica_min, input$cica_max), 0, 1), clamp(max(input$cica_min, input$cica_max), 0, 1))
  })
  ylim_wue <- reactive({
    if (!isTRUE(input$fix_axes)) return(NULL)
    c(min(input$wue_min, input$wue_max), max(input$wue_min, input$wue_max))
  })
  
  
  # -------------------------
  # Treatment presets
  # -------------------------
  observeEvent(input$preset, {
    if (input$preset == "Ambient (A)") {
      updateSliderInput(session, "Tair_C", value = 25)
      updateSliderInput(session, "REW", value = 1)
    } else if (input$preset == "Heat (H)") {
      updateSliderInput(session, "Tair_C", value = 30)
      updateSliderInput(session, "REW", value = 1)
    } else if (input$preset == "Drought (D)") {
      updateSliderInput(session, "Tair_C", value = 25)
      updateSliderInput(session, "REW", value = 0.4)
    } else if (input$preset == "Heat + Drought (H+D)") {
      updateSliderInput(session, "Tair_C", value = 30)
      updateSliderInput(session, "REW", value = 0.4)
    }
  })
  
  lambda_eff <- reactive({
    
    ca_molmol <- input$ca_ppm * 1e-6
    
    # baseline λ from inference or slider
    if (!isTRUE(input$infer_lambda)) {
      lambda_base <- pmax(input$lambda_in, 1e-12)
    } else {
      Dref <- pmax(input$VPD_ref / input$patm_kPa, 1e-12)
      lambda_base <- lambda_from_ci_ca_eq20(
        ci_ca = input$ci_ca_ref,
        D_molmol = Dref,
        ca_molmol = ca_molmol,
        a = input$a
      )
      lambda_base <- pmax(lambda_base, 1e-12)
    }
    
    # choose wet reference consistently
    lambda_wet <- if (isTRUE(input$infer_lambda)) {
      lambda_base
    } else {
      pmax(input$lambda_wet, 1e-12)
    }
    
    # apply drought scaling
    lam <- lambda_from_REW(
      REW = input$REW,
      lambda_wet = lambda_wet,
      form = input$drought_fun,
      k = input$k_rew,
      beta = input$beta_rew
    )
    
    pmax(lam, 1e-12)
  })
  
  model_df <- reactive({
    vmin <- min(input$VPDmin, input$VPDmax)
    vmax <- max(input$VPDmin, input$VPDmax)
    
    VPD_kPa <- seq(vmin, vmax, length.out = input$npts)
    VPD_kPa <- pmax(VPD_kPa, 0)
    
    D_molmol <- pmax(VPD_kPa / input$patm_kPa, 1e-12)
    
    # Convert ppm -> mol/mol
    ca_molmol <- input$ca_ppm * 1e-6
    alpha2_molmol <- input$alpha2_ppm * 1e-6
    cp_molmol <- input$GammaStar * 1e-6
    
    # alpha2' used in Eq. 16-like term (paper defines alpha2' = alpha2 + ca)
    alpha2_prime_molmol <- alpha2_molmol + ca_molmol
    
    # Convert alpha1 µmol -> mol
    # alpha1_mol <- input$alpha1_umol * 1e-6
    
    # Now alpha1 with heat effect
    alpha1_umol_eff <- if (isTRUE(input$use_T_for_alpha1)) {
      vcmax_T_table1(input$Tair_C, input$Vcmax25, input$m1, input$m2)
    } else {
      input$alpha1_umol
    }
    alpha1_mol <- alpha1_umol_eff * 1e-6
    
    lam <- lambda_eff()
    
    gco2_mol <- g_opt(
      D_molmol = D_molmol, ca_molmol = ca_molmol,
      a = input$a, lambda = lam,
      alpha1_mol = alpha1_mol, alpha2_molmol = alpha2_molmol, s = input$s
    )
    
    open <- gco2_mol > 0
    
    gh2o_mol <- input$a * gco2_mol
    
    fc_mol <- fc_from_g(
      gco2_mol,
      ca_molmol = ca_molmol, alpha1_mol = alpha1_mol,
      alpha2_molmol = alpha2_molmol, s = input$s
    )
    
    fe_mol <- fe_from_g(gco2_mol, D_molmol = D_molmol, a = input$a)
    
    # Conversions for plotting
    fc_umol <- fc_mol * 1e6   # µmol CO2 m^-2 s^-1
    fe_mmol <- fe_mol * 1e3   # mmol H2O m^-2 s^-1
    
    # ci/ca from Eq. 20 in open-stomata regime; otherwise ci/ca -> 1 (physical limit)
    ci_ca_eq20 <- ci_ca_opt(D_molmol = D_molmol, ca_molmol = ca_molmol, a = input$a, lambda = lam)
    ci_ca <- ifelse(open, ci_ca_eq20, 1.0)
    
    # WUE: flux ratio and Eq. 21
    wue_flux <- ifelse(fe_mol > 0, fc_mol / fe_mol, NA_real_)
    wue_eq21 <- wue_closed_eq21(D_molmol = D_molmol, ca_molmol = ca_molmol, a = input$a, lambda = lam)
    
    # AFM diagnostic derivative: use only open regime to avoid clamp artifacts
    fe_mmol_open <- ifelse(open, fe_mmol, NA_real_)
    dfe_dVPD <- dydx_fd(VPD_kPa, fe_mmol_open)
    
    # numeric peak in selected range (open regime)
    fe_for_peak <- ifelse(open, fe_mmol, NA_real_)
    i_peak <- which.max(fe_for_peak)
    VPD_peak <- if (length(i_peak) == 0 || is.infinite(i_peak) || is.na(fe_for_peak[i_peak])) NA_real_ else VPD_kPa[i_peak]
    
    # TRUE marginal WUE (Nakad/Katul): ∂fc/∂fe|D = (∂fc/∂g)/(∂fe/∂g) = λ
    dfc_dg_val <- dfc_dg(gco2_mol, ca_molmol, alpha1_mol, alpha2_molmol, input$s)
    mwue_true <- ifelse(open, dfc_dg_val / dfe_dg(D_molmol, input$a), NA_real_)
    mwue_resid <- mwue_true - lam
    
    data.frame(
      VPD_kPa = VPD_kPa,
      D_molmol = D_molmol,
      open = open,
      gco2_mol = gco2_mol,
      gh2o_mol = gh2o_mol,
      fc_umol = fc_umol,
      fe_mmol = fe_mmol,
      ci_ca = ci_ca,
      wue = wue_flux,
      wue_closed = wue_eq21,
      dfe_dVPD = dfe_dVPD,
      VPD_peak = VPD_peak,
      mwue_true = mwue_true,
      mwue_resid = mwue_resid,
      alpha2_prime_molmol = alpha2_prime_molmol,
      cp_molmol = cp_molmol
    )
  })
  
  Dafm_kPa <- reactive({
    ca_molmol <- input$ca_ppm * 1e-6
    lam <- lambda_eff()
    D_afm_molmol(ca_molmol = ca_molmol, a = input$a, lambda = lam) * input$patm_kPa
  })
  
  Dafm_D <- reactive({
    ca_molmol <- input$ca_ppm * 1e-6
    lam <- lambda_eff()
    D_afm_molmol(ca_molmol = ca_molmol, a = input$a, lambda = lam)
  })
  
  # -------------------------
  # Plots: Curves
  # -------------------------
  output$p_gco2 <- renderPlot({
    df <- model_df()
    xafm <- Dafm_kPa()
    xmin <- min(df$VPD_kPa); xmax <- max(df$VPD_kPa)
    
    ggplot(df, aes(VPD_kPa, gco2_mol)) +
      geom_line() +
      vline_if_in_range(xafm, xmin, xmax, "dashed") +
      coord_cartesian(xlim = c(xmin, xmax), ylim = ylim_gco2()) +
      labs(
        x = "VPD (kPa)",
        y = expression(paste(g[CO[2]], " (mol CO"[2], " m"^{-2}, " s"^{-1}, ")")),
        title = "Optimal stomatal conductance to CO₂ vs VPD (dashed: theory Dafm)"
      )
  })
  
  output$p_gh2o <- renderPlot({
    df <- model_df()
    xafm <- Dafm_kPa()
    xmin <- min(df$VPD_kPa); xmax <- max(df$VPD_kPa)
    
    ggplot(df, aes(VPD_kPa, gh2o_mol)) +
      geom_line() +
      vline_if_in_range(xafm, xmin, xmax, "dashed") +
      coord_cartesian(xlim = c(xmin, xmax), ylim = ylim_gh2o()) +
      labs(
        x = "VPD (kPa)",
        y = expression(paste(g[H[2]*O], " (mol H"[2], "O m"^{-2}, " s"^{-1}, ")")),
        title = "Equivalent stomatal conductance to H₂O vs VPD (gH₂O = a·gCO₂)"
      )
  })
  
  output$p_fc <- renderPlot({
    df <- model_df()
    xafm <- Dafm_kPa()
    xmin <- min(df$VPD_kPa); xmax <- max(df$VPD_kPa)
    
    ggplot(df, aes(VPD_kPa, fc_umol)) +
      geom_line() +
      vline_if_in_range(xafm, xmin, xmax, "dashed") +
      coord_cartesian(xlim = c(xmin, xmax), ylim = ylim_fc()) +
      labs(
        x = "VPD (kPa)",
        y = expression(paste(f[c], " (", mu, "mol CO"[2], " m"^{-2}, " s"^{-1}, ")")),
        title = "Photosynthesis vs VPD (with optimal g)"
      )
  })
  
  output$p_fe <- renderPlot({
    df <- model_df()
    xafm <- Dafm_kPa()
    xpeak <- df$VPD_peak[1]
    xmin <- min(df$VPD_kPa); xmax <- max(df$VPD_kPa)
    
    ggplot(df, aes(VPD_kPa, fe_mmol)) +
      geom_line() +
      vline_if_in_range(xafm, xmin, xmax, "dashed") +
      vline_if_in_range(xpeak, xmin, xmax, "dotted") +
      coord_cartesian(xlim = c(xmin, xmax), ylim = ylim_fe()) +
      labs(
        x = "VPD (kPa)",
        y = expression(paste(f[e], " (mmol H"[2], "O m"^{-2}, " s"^{-1}, ")")),
        title = "Transpiration proxy vs VPD (dashed: theory Dafm; dotted: numeric peak)"
      )
  })
  
  output$p_ci <- renderPlot({
    df <- model_df()
    xafm <- Dafm_kPa()
    xmin <- min(df$VPD_kPa); xmax <- max(df$VPD_kPa)
    
    ggplot(df, aes(VPD_kPa, ci_ca)) +
      geom_line() +
      geom_hline(yintercept = 1, linetype = "dotted") +
      vline_if_in_range(xafm, xmin, xmax, "dashed") +
      coord_cartesian(xlim = c(xmin, xmax), ylim = ylim_cica()) +
      labs(
        x = "VPD (kPa)",
        y = expression(paste(c[i] / c[a], " (-)")),
        title = "ci/ca vs VPD (Eq. 20 in open regime; ci/ca→1 when g=0)"
      )
  })
  
  output$p_wue <- renderPlot({
    df <- model_df()
    xafm <- Dafm_kPa()
    xmin <- min(df$VPD_kPa); xmax <- max(df$VPD_kPa)
    lam <- lambda_eff()
    
    ggplot(df, aes(VPD_kPa)) +
      geom_line(aes(y = wue), size = 1) +
      geom_line(aes(y = wue_closed), linetype = "dotted", size = 1) +
      vline_if_in_range(xafm, xmin, xmax, "dashed") +
      
      # Horizontal line showing λ (true marginal WUE)
      geom_hline(yintercept = lam, linetype = "dotdash") +
      
      # Text annotation with λ value
      annotate("text",
               x = xmin + 0.05*(xmax - xmin),
               y = lam,
               label = paste0("mWUE = λ = ", signif(lam, 4)),
               hjust = 0,
               vjust = -0.5,
               size = 5) +
      
      coord_cartesian(xlim = c(xmin, xmax), ylim = ylim_wue()) +
      labs(
        x = "VPD (kPa)",
        y = expression(paste("WUE (mol CO"[2], " / mol H"[2], "O)")),
        title = "WUE vs VPD (solid: fc/fe; dotted: Eq. 21; dotdash: λ = marginal WUE)"
      ) +
      theme(
        axis.title = element_text(size = 14),
        axis.text  = element_text(size = 12),
        plot.title = element_text(size = 14)
      )
  })
  
  # -------------------------
  # AFM diagnostics
  # -------------------------
  output$p_dfe <- renderPlot({
    df <- model_df()
    xafm <- Dafm_kPa()
    xpeak <- df$VPD_peak[1]
    xmin <- min(df$VPD_kPa); xmax <- max(df$VPD_kPa)
    
    ggplot(df, aes(VPD_kPa, dfe_dVPD)) +
      geom_line() +
      geom_hline(yintercept = 0, linetype = "dashed") +
      vline_if_in_range(xafm, xmin, xmax, "dashed") +
      vline_if_in_range(xpeak, xmin, xmax, "dotted") +
      coord_cartesian(xlim = c(xmin, xmax)) +
      labs(
        x = "VPD (kPa)",
        y = expression(paste("d f"[e], "/ d(VPD) (mmol m"^{-2}, " s"^{-1}, " kPa"^{-1}, ")")),
        title = "AFM diagnostic (open-regime derivative; dashed: theory Dafm; dotted: numeric peak)"
      )
  })
  
  output$txt_afm <- renderPrint({
    ca_molmol <- input$ca_ppm * 1e-6
    alpha2_molmol <- input$alpha2_ppm * 1e-6
    cp_molmol <- input$GammaStar * 1e-6
    alpha2_prime_molmol <- alpha2_molmol + ca_molmol
    lam <- lambda_eff()
    
    Dafm_m <- Dafm_D()
    Dafm_k <- Dafm_m * input$patm_kPa
    
    cafm <- Ca_afm_coeff(alpha2_prime_molmol = alpha2_prime_molmol, ca_molmol = ca_molmol, cp_molmol = cp_molmol)
    
    df <- model_df()
    VPD_peak <- df$VPD_peak[1]
    fe_peak <- suppressWarnings(max(ifelse(df$open, df$fe_mmol, NA_real_), na.rm = TRUE))
    if (!is.finite(fe_peak)) fe_peak <- NA_real_
    
    idx <- which(!is.na(df$dfe_dVPD) & df$dfe_dVPD <= 0)
    VPD_cross <- if (length(idx) > 0) df$VPD_kPa[min(idx)] else NA_real_
    D_cross <- if (!is.na(VPD_cross)) VPD_cross / input$patm_kPa else NA_real_
    
    frac_open <- mean(df$open)
    
    cat("Mode:\n")
    cat(if (isTRUE(input$infer_lambda)) "  Diagnostic: λ inferred from Eq. 20 inversion\n" else "  Forward: λ set by slider\n")
    cat(sprintf("  λ used = %.6g\n\n", lam))
    
    if (isTRUE(input$infer_lambda)) {
      cat(sprintf("  ci/ca(ref) = %.3f at VPDref = %.2f kPa\n", input$ci_ca_ref, input$VPD_ref))
      cat(sprintf("  Dref = VPDref/Patm = %.6g (mol/mol)\n\n", input$VPD_ref / input$patm_kPa))
    }
    
    cat("AFM lines and peak check:\n")
    cat(sprintf("  Theory Dafm (D units) = %.6g (mol/mol)\n", Dafm_m))
    cat(sprintf("  Theory Dafm (kPa)     = %.3f kPa\n", Dafm_k))
    cat(sprintf("  Numeric peak VPD      = %s kPa (max fe(open) = %s mmol m^-2 s^-1)\n",
                ifelse(is.na(VPD_peak), "NA", sprintf("%.3f", VPD_peak)),
                ifelse(is.na(fe_peak), "NA", sprintf("%.3f", fe_peak))))
    
    if (is.na(VPD_cross)) {
      cat("  First d(fe)/dVPD <= 0  = not found in selected range\n")
    } else {
      cat(sprintf("  First d(fe)/dVPD <= 0  = %.3f kPa (D = %.6g mol/mol)\n", VPD_cross, D_cross))
    }
    
    cat("\nOpen-stomata regime check:\n")
    cat(sprintf("  Fraction of VPD grid with g>0 = %.2f\n", frac_open))
    if (frac_open < 0.25) {
      cat("  WARNING: Most of the chosen VPD range is in the g=0 regime; theory comparisons may look odd.\n")
    }
    
    cat("\nWhy theory Dafm may not equal numeric peak:\n")
    cat("- Dafm is a simplified analytical onset estimate; fe uses full g(D) with truncation g>=0.\n")
    cat("- The peak is computed within a finite VPD window.\n")
    cat("- Therefore we show both (dashed: theory; dotted: numeric peak).\n\n")
    
    cat(sprintf("Eq.16-like coefficient Ca^afm = %.6g  (paper discussion uses Ca^afm < 0)\n", cafm))
  })
  
  output$txt_mwue <- renderPrint({
    df <- model_df()
    lam <- lambda_eff()
    
    resid <- df$mwue_resid[df$open]
    resid <- resid[is.finite(resid)]
    
    cat("Marginal WUE (Nakad/Katul definition):\n")
    cat("  mWUE ≡ (∂f_c/∂f_e)|D = (∂f_c/∂g)/(∂f_e/∂g) = λ\n\n")
    cat(sprintf("  λ used = %.6g\n", lam))
    
    if (length(resid) == 0) {
      cat("\nConsistency check:\n")
      cat("  No open-stomata points (g>0) in the selected VPD range.\n")
      cat("  Try reducing VPDmax or λ, or increasing alpha1.\n")
      return(invisible(NULL))
    }
    
    cat("\nConsistency check (open-stomata regime):\n")
    cat(sprintf("  max|mWUE - λ| = %.3e\n", max(abs(resid))))
    cat(sprintf("  mean(mWUE - λ) = %.3e\n", mean(resid)))
    cat(sprintf("  sd(mWUE - λ)   = %.3e\n", sd(resid)))
    
    cat("\nInterpretation:\n")
    cat("  This is a local (partial) marginal tradeoff with respect to g at fixed D.\n")
    cat("  It is the quantity referred to as marginal WUE in Katul et al. (2012) and Nakad et al. (2023).\n")
  })
  
  # -------------------------
  # A–Ci (Farquhar) page
  # -------------------------
  output$p_aci <- renderPlot({
    # Ci range for curve (0..~1.2*Ca)
    Ci_ppm <- seq(0, max(10, 1.2 * input$ca_ppm), length.out = 400)
    Ci_ppm <- pmax(Ci_ppm, 1e-6)
    
    aci <- farquhar_aci(
      Ci_ppm = Ci_ppm,
      Vcmax = input$Vcmax,
      J = input$J,
      Rd = input$Rd,
      GammaStar_ppm = input$GammaStar,
      Kc_ppm = input$Kc,
      Ko_ppm = input$Ko,
      O2_ppm = input$O2
    )
    
    # Operating Ci predicted by Nakad-optimality at chosen VPD
    VPDk <- clamp(input$VPD_for_aci, min(input$VPDmin, input$VPDmax), max(input$VPDmin, input$VPDmax))
    Dk <- pmax(VPDk / input$patm_kPa, 1e-12)
    ca_molmol <- input$ca_ppm * 1e-6
    lam <- lambda_eff()
    ci_ca_k <- ci_ca_opt(D_molmol = Dk, ca_molmol = ca_molmol, a = input$a, lambda = lam)
    ci_ca_k <- clamp(ci_ca_k, 0, 1)
    Ci_oper_ppm <- ci_ca_k * input$ca_ppm
    
    dd <- rbind(
      data.frame(Ci_ppm = aci$Ci_ppm, A = aci$A, type = "A = min(Ac,Aj) - Rd"),
      data.frame(Ci_ppm = aci$Ci_ppm, A = aci$Ac, type = "Ac (Rubisco-limited)"),
      data.frame(Ci_ppm = aci$Ci_ppm, A = aci$Aj, type = "Aj (RuBP-limited)")
    )
    
    ggplot(dd, aes(Ci_ppm, A, color = type, linetype = type)) +
      geom_line(size = 1) +
      scale_color_manual(values = c(
        "A = min(Ac,Aj) - Rd" = "black",
        "Ac (Rubisco-limited)" = "grey20",
        "Aj (RuBP-limited)" = "red"
      )) +
      geom_vline(xintercept = Ci_oper_ppm, linetype = "dashed") +
      labs(
        x = expression(paste(C[i], " (ppm)")),
        y = expression(paste("A (", mu, "mol m"^{-2}, " s"^{-1}, ")")),
        title = sprintf("A–Ci (Farquhar) with operating Ci marker (VPD = %.2f kPa)", VPDk)
      ) +
      theme(
        legend.position = "bottom",
        legend.text = element_text(size = 14),
        legend.title = element_blank(),
        axis.title = element_text(size = 14),
        axis.text  = element_text(size = 12),
        plot.title = element_text(size = 14)
      )
  })
  
  # -------------------------
  # WUE geometry page
  # -------------------------
  output$p_wue_geom <- renderPlot({
    df <- model_df()
    lam <- lambda_eff()
    
    df_open <- df[df$open & is.finite(df$fe_mmol) & is.finite(df$fc_umol) & df$fe_mmol > 0, , drop=FALSE]
    if (nrow(df_open) < 10) { plot.new(); text(0.5,0.5,"Not enough open points"); return() }
    
    # split at fe peak (avoids looping comb)
    i_peak <- which.max(df_open$fe_mmol)
    if (is.na(i_peak) || i_peak < 3) i_peak <- floor(nrow(df_open)/2)
    df1 <- df_open[1:i_peak, , drop=FALSE]
    df2 <- df_open[i_peak:nrow(df_open), , drop=FALSE]
    
    # choose tangent point on rising branch (pre-AFM)
    idx_mid <- max(2, floor(nrow(df1)/2))
    fe0 <- df1$fe_mmol[idx_mid]
    fc0 <- df1$fc_umol[idx_mid]
    
    wue_avg <- fc0 / fe0
    
    fe_range <- range(df_open$fe_mmol, na.rm=TRUE)
    fe_line <- seq(fe_range[1], fe_range[2], length.out = 100)
    
    lambda_plot <- lam * 1000  # mol/mol -> µmol/mmol
    fc_tangent <- fc0 + lambda_plot * (fe_line - fe0)
    fc_avgline <- wue_avg * fe_line
    
    ggplot() +
      geom_line(data = df1, aes(fe_mmol, fc_umol), linewidth = 1.2) +
      geom_line(data = df2, aes(fe_mmol, fc_umol), linewidth = 1.2, linetype = "dotted") +
      geom_point(data = data.frame(fe0=fe0, fc0=fc0), aes(fe0, fc0), size = 3) +
      geom_line(data = data.frame(fe_line=fe_line, fc_line=fc_avgline),
                aes(fe_line, fc_line), inherit.aes=FALSE, linetype="dashed", linewidth=1) +
      geom_line(data = data.frame(fe_line=fe_line, fc_line=fc_tangent),
                aes(fe_line, fc_line), inherit.aes=FALSE, linewidth=1.2) +
      labs(
        x = expression(paste(f[e], " (mmol H"[2], "O m"^{-2}, " s"^{-1}, ")")),
        y = expression(paste(f[c], " (", mu, "mol CO"[2], " m"^{-2}, " s"^{-1}, ")")),
        title = "WUE geometry: average vs marginal (solid: pre-AFM, dotted: post-AFM)"
      )
  })
}

shinyApp(ui, server)