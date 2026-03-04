# app.R
# Nakad et al. (2023): Stomatal optimality, VPD response, and AFM (exploratory tool)
#
# This version is integrated + internally consistent with the *simplified* Nakad framework:
#   - Eq. 9 (Fick: fc = g(ca-ci), fe ≈ a g D), Eq. 10 (linearized demand form used here),
#     Eq. 12 (optimal g), Eq. 20 (optimal ci/ca), Eq. 21 (WUE), and simplified AFM Dafm expression.
#   - Diagnostic mode: infer λ by direct inversion of Eq. 20 at VPDref and target ci/ca.
#   - Adds H2O conductance: g_H2O = a * g_CO2.
#   - Shows both: dashed = theory Dafm; dotted = numeric peak of fe within chosen VPD range.
#   - Adds optional fixed y-axis limits for all curves (user-controlled).
#
# Unit conventions (internally consistent):
#   - ca, alpha2, cp entered in ppm; converted to mol/mol using *1e-6
#   - alpha1 entered in µmol m^-2 s^-1; converted to mol m^-2 s^-1 using *1e-6
#   - VPD entered in kPa; converted to D (dimensionless mol/mol) via D = VPD/Patm
#   - g_CO2 in mol CO2 m^-2 s^-1; g_H2O in mol H2O m^-2 s^-1 (equivalent conductance)
#   - fc plotted in µmol CO2 m^-2 s^-1; fe plotted in mmol H2O m^-2 s^-1
#   - WUE in mol CO2 / mol H2O

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

# Eq. 20: optimal ci/ca (dimensionless)
ci_ca_opt <- function(D_molmol, ca_molmol, a, lambda) {
  D_molmol <- pmax(D_molmol, 1e-12)
  lambda <- pmax(lambda, 1e-12)
  1 - sqrt(a * lambda / ca_molmol) * sqrt(D_molmol)
}

# Invert Eq. 20 exactly to infer lambda from (ci/ca at Dref, ca):
# ci/ca = 1 - sqrt(a*lambda/ca)*sqrt(D)
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

# Eq. 16-like sign term (used in discussion; inputs in mol/mol)
Ca_afm_coeff <- function(alpha2_prime_molmol, ca_molmol, cp_molmol) {
  - (alpha2_prime_molmol - 2 * ca_molmol + 2 * cp_molmol) / (alpha2_prime_molmol^2)
}

# numerical derivative d(fe)/d(VPD) (finite differences), VPD in kPa, fe in mmol m^-2 s^-1
dfe_dVPD_numeric <- function(VPD_kPa, fe_mmol) {
  dX <- diff(VPD_kPa)
  dY <- diff(fe_mmol)
  c(NA, dY / dX)
}

# -------------------------
# UI
# -------------------------

ui <- fluidPage(
  titlePanel("Nakad et al. (2023): Stomatal optimality, VPD response, and AFM (exploratory tool)"),
  
  sidebarLayout(
    sidebarPanel(
      h4("Inputs"),
      
      sliderInput("ca_ppm", "Atmospheric CO₂, ca (ppm)", min = 250, max = 1000, value = 420, step = 10),
      
      checkboxInput("infer_lambda", "Infer λ from ci/ca (diagnostic mode; Eq. 20 inversion)", value = TRUE),
      
      conditionalPanel(
        condition = "input.infer_lambda == false",
        sliderInput("lambda_in", "λ (input; dimensionless in this simplified form)",
                    min = 1e-5, max = 0.5, value = 0.03, step = 0.001)
      ),
      
      conditionalPanel(
        condition = "input.infer_lambda == true",
        sliderInput("ci_ca_ref", HTML("Reference c<sub>i</sub>/c<sub>a</sub> at VPD<sub>ref</sub>"),
                    min = 0.2, max = 0.95, value = 0.7, step = 0.01),
        sliderInput("VPD_ref", HTML("VPD<sub>ref</sub> (kPa) used to infer &lambda;"),
                    min = 0.1, max = 5, value = 1.0, step = 0.1)
      ),
      
      sliderInput(
        "alpha1_umol",
        label = HTML("&alpha;<sub>1</sub> (&mu;mol CO<sub>2</sub> m<sup>-2</sup> s<sup>-1</sup>; Vcmax-like)"),
        min = 10, max = 200, value = 80, step = 1
      ),
      
      sliderInput(
        "alpha2_ppm",
        label = HTML("&alpha;<sub>2</sub> (ppm; e.g., K<sub>c</sub>(1+O/K<sub>o</sub>) under RuBisCO limitation)"),
        min = 50, max = 1200, value = 550, step = 10
      ),
      
      sliderInput("s", "s (ci/ca linearization parameter)", min = 0.3, max = 0.9, value = 0.7, step = 0.01),
      sliderInput("a", "a (H₂O/CO₂ diffusivity ratio)", min = 1.2, max = 2.0, value = 1.6, step = 0.01),
      
      sliderInput(
        "cp_ppm",
        label = HTML("c<sub>p</sub> (CO<sub>2</sub> compensation point, ppm)"),
        min = 0, max = 100, value = 40, step = 1
      ),
      
      sliderInput("patm_kPa", "Atmospheric pressure, Patm (kPa)", min = 80, max = 110, value = 101.3, step = 0.1),
      
      hr(),
      h4("VPD range"),
      sliderInput("VPDmin", "VPD min (kPa)", min = 0, max = 5, value = 0.1, step = 0.05),
      sliderInput("VPDmax", "VPD max (kPa)", min = 0.2, max = 10, value = 5, step = 0.1),
      numericInput("npts", "Number of points", value = 500, min = 100, max = 2000, step = 50),
      
      hr(),
      h4("Axis limits"),
      checkboxInput("fix_axes", "Fix y-axis limits (for comparisons)", value = FALSE),
      
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
          numericInput("fe_max", "max (mmol H2O m^-2 s^-1)", value = 20, min = 0),
          
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
      h4("Notes"),
      helpText(
        "We convert VPD (kPa) to D = VPD/Patm (dimensionless mol/mol) to keep Eq. 9 unit-consistent.",
        "Transpiration proxy: fe ≈ a g D.",
        "On fe plot: dashed = theory Dafm; dotted = numeric peak of fe in chosen range."
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
                 verbatimTextOutput("txt_afm")
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
  
  lambda_eff <- reactive({
    ca_molmol <- input$ca_ppm * 1e-6
    
    if (!isTRUE(input$infer_lambda)) {
      return(pmax(input$lambda_in, 1e-12))
    }
    
    Dref <- pmax(input$VPD_ref / input$patm_kPa, 1e-12)
    lam <- lambda_from_ci_ca_eq20(
      ci_ca = input$ci_ca_ref,
      D_molmol = Dref,
      ca_molmol = ca_molmol,
      a = input$a
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
    cp_molmol <- input$cp_ppm * 1e-6
    
    # alpha2' used in Eq. 16-like term (paper defines alpha2' = alpha2 + ca)
    alpha2_prime_molmol <- alpha2_molmol + ca_molmol
    
    # Convert alpha1 µmol -> mol
    alpha1_mol <- input$alpha1_umol * 1e-6
    
    lam <- lambda_eff()
    
    gco2_mol <- g_opt(
      D_molmol = D_molmol, ca_molmol = ca_molmol,
      a = input$a, lambda = lam,
      alpha1_mol = alpha1_mol, alpha2_molmol = alpha2_molmol, s = input$s
    )
    
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
    
    # ci/ca from Eq. 20
    ci_ca_eq20 <- ci_ca_opt(D_molmol = D_molmol, ca_molmol = ca_molmol, a = input$a, lambda = lam)
    
    # WUE: flux ratio and Eq. 21
    wue_flux <- ifelse(fe_mol > 0, fc_mol / fe_mol, NA_real_)
    wue_eq21 <- wue_closed_eq21(D_molmol = D_molmol, ca_molmol = ca_molmol, a = input$a, lambda = lam)
    
    dfe_dVPD <- dfe_dVPD_numeric(VPD_kPa, fe_mmol)
    
    # numeric peak in selected range
    i_peak <- which.max(fe_mmol)
    VPD_peak <- VPD_kPa[i_peak]
    
    data.frame(
      VPD_kPa = VPD_kPa,
      D_molmol = D_molmol,
      gco2_mol = gco2_mol,
      gh2o_mol = gh2o_mol,
      fc_umol = fc_umol,
      fe_mmol = fe_mmol,
      ci_ca = ci_ca_eq20,
      wue = wue_flux,
      wue_closed = wue_eq21,
      dfe_dVPD = dfe_dVPD,
      VPD_peak = VPD_peak,
      alpha2_prime_molmol = alpha2_prime_molmol,
      cp_molmol = cp_molmol
    )
  })
  
  Dafm_kPa <- reactive({
    ca_molmol <- input$ca_ppm * 1e-6
    lam <- lambda_eff()
    Dafm_m <- D_afm_molmol(ca_molmol = ca_molmol, a = input$a, lambda = lam)
    Dafm_m * input$patm_kPa
  })
  
  # -------------------------
  # Plots
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
      geom_hline(yintercept = 0, linetype = "dotted") +
      vline_if_in_range(xafm, xmin, xmax, "dashed") +
      coord_cartesian(xlim = c(xmin, xmax), ylim = ylim_cica()) +
      labs(
        x = "VPD (kPa)",
        y = expression(paste(c[i] / c[a], " (-)")),
        title = "Optimal ci/ca vs VPD (Eq. 20)"
      )
  })
  
  output$p_wue <- renderPlot({
    df <- model_df()
    xafm <- Dafm_kPa()
    xmin <- min(df$VPD_kPa); xmax <- max(df$VPD_kPa)
    
    ggplot(df, aes(VPD_kPa)) +
      geom_line(aes(y = wue)) +
      geom_line(aes(y = wue_closed), linetype = "dotted") +
      vline_if_in_range(xafm, xmin, xmax, "dashed") +
      coord_cartesian(xlim = c(xmin, xmax), ylim = ylim_wue()) +
      labs(
        x = "VPD (kPa)",
        y = expression(paste("WUE (mol CO"[2], " / mol H"[2], "O)")),
        title = "WUE vs VPD (solid: fc/fe; dotted: Eq. 21)"
      )
  })
  
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
        title = "AFM diagnostic (dashed: theory Dafm; dotted: numeric peak)"
      )
  })
  
  output$txt_afm <- renderPrint({
    ca_molmol <- input$ca_ppm * 1e-6
    alpha2_molmol <- input$alpha2_ppm * 1e-6
    cp_molmol <- input$cp_ppm * 1e-6
    alpha2_prime_molmol <- alpha2_molmol + ca_molmol
    lam <- lambda_eff()
    
    Dafm_m <- D_afm_molmol(ca_molmol = ca_molmol, a = input$a, lambda = lam)
    Dafm_k <- Dafm_m * input$patm_kPa
    
    cafm <- Ca_afm_coeff(alpha2_prime_molmol = alpha2_prime_molmol, ca_molmol = ca_molmol, cp_molmol = cp_molmol)
    
    df <- model_df()
    VPD_peak <- df$VPD_peak[1]
    fe_peak <- max(df$fe_mmol, na.rm = TRUE)
    
    idx <- which(!is.na(df$dfe_dVPD) & df$dfe_dVPD <= 0)
    VPD_cross <- if (length(idx) > 0) df$VPD_kPa[min(idx)] else NA_real_
    
    cat("Mode:\n")
    cat(if (isTRUE(input$infer_lambda)) "  Diagnostic: λ inferred from Eq. 20 inversion\n" else "  Forward: λ set by slider\n")
    cat(sprintf("  λ used = %.6g\n\n", lam))
    
    if (isTRUE(input$infer_lambda)) {
      cat(sprintf("  ci/ca(ref) = %.3f at VPDref = %.2f kPa\n", input$ci_ca_ref, input$VPD_ref))
      cat(sprintf("  (Dref = VPDref/Patm = %.6g)\n\n", input$VPD_ref / input$patm_kPa))
    }
    
    cat("AFM lines and peak check:\n")
    cat(sprintf("  Theory Dafm (D-units) = %.6g (mol/mol)\n", Dafm_m))
    cat(sprintf("  Theory Dafm (kPa)     = %.3f kPa\n", Dafm_k))
    cat(sprintf("  Numeric peak VPD      = %.3f kPa (max fe = %.3f mmol m^-2 s^-1)\n", VPD_peak, fe_peak))
    
    if (is.na(VPD_cross)) {
      cat("  First d(fe)/dVPD <= 0  = not found in selected range\n\n")
    } else {
      cat(sprintf("  First d(fe)/dVPD <= 0  = %.3f kPa\n\n", VPD_cross))
    }
    
    cat("Why theory Dafm may not equal numeric peak:\n")
    cat("- Dafm is a simplified analytical onset estimate; plotted fe uses full g(D) expression plus truncation g>=0.\n")
    cat("- The peak is computed within a finite VPD window.\n")
    cat("- Therefore we show both (dashed: theory; dotted: numeric peak).\n\n")
    
    cat(sprintf("Eq.16-like coefficient Ca^afm = %.6g  (paper discussion uses Ca^afm < 0)\n", cafm))
  })
}

shinyApp(ui, server)