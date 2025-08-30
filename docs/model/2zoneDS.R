library(shiny)
library(ggplot2)
library(dplyr)
library(shinydashboard)

# 定義 UI
ui <- dashboardPage(
  dashboardHeader(title = "二暴露區模式暴露評估系統"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("輸入參數", tabName = "input", icon = icon("sliders")),
      menuItem("計算結果", tabName = "result", icon = icon("chart-line")),
      menuItem("風險分級", tabName = "risk", icon = icon("exclamation-triangle")),
      menuItem("圖表展示", tabName = "plot", icon = icon("chart-area"))
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName = "input",
              fluidRow(
                box(title = "基本參數", width = 6,
                    numericInput("G", "污染物逸散率 G (mg/min)", value = 1000),
                    numericInput("Q", "房間通風率 Q (m³/min)", value = 20),
                    numericInput("beta", "近遠場交換率 β (m³/min)", value = 6.7),
                    numericInput("V_N", "近場體積 V_N (m³)", value = 0.93),
                    numericInput("V_F", "遠場體積 V_F (m³)", value = 199),
                    numericInput("t_work", "作業時間 (分鐘)", value = 30),
                    numericInput("t_decay_start", "衰減開始時間 (分鐘)", value = 30),
                    numericInput("t_total", "總模擬時間 (分鐘)", value = 60)
                ),
                box(title = "暴露限值設定", width = 6,
                    selectInput("unit", "單位", choices = c("mg/m³", "ppm")),
                    numericInput("TWA", "TWA (8hr) 限值", value = 50),
                    numericInput("STEL", "STEL (15min) 限值", value = 100),
                    numericInput("IDLH", "IDLH 限值", value = 500),
                    numericInput("OEL", "OEL 限值", value = 50)
                )
              ),
              fluidRow(
                box(title = "功能選擇", width = 12,
                    selectInput("func", "選擇計算功能", 
                                choices = c("時間濃度解", "穩態濃度", "濃度衰減")),
                    actionButton("calc", "開始計算", icon = icon("calculator"))
                )
              )
      ),
      tabItem(tabName = "result",
              fluidRow(
                box(title = "計算結果", width = 12,
                    verbatimTextOutput("result_text")
                )
              )
      ),
      tabItem(tabName = "risk",
              fluidRow(
                box(title = "暴露風險分級", width = 12,
                    tableOutput("risk_table")
                )
              )
      ),
      tabItem(tabName = "plot",
              fluidRow(
                box(title = "濃度趨勢圖", width = 12,
                    plotOutput("conc_plot")
                )
              )
      )
    )
  )
)

# 定義 Server
server <- function(input, output) {
  
  # 計算 λ1, λ2
  calc_lambda <- function(beta, Q, V_N, V_F) {
    A <- beta * V_F + V_N * (beta + Q)
    B <- V_N * V_F
    C <- 4 * beta * Q / B
    D <- (A / B)^2 - C
    lambda1 <- 0.5 * (-A/B + sqrt(D))
    lambda2 <- 0.5 * (-A/B - sqrt(D))
    return(list(lambda1 = lambda1, lambda2 = lambda2))
  }
  
  # 公式6-3: 時間濃度解 - 近場
  calc_C_N_time <- function(t, G, beta, Q, V_N, V_F, lambda1, lambda2) {
    term1 <- G / (beta * Q / (beta + Q))
    term2 <- G * (beta * Q + lambda2 * V_N * (beta + Q)) / (beta * Q * V_N * (lambda1 - lambda2)) * exp(lambda1 * t)
    term3 <- G * (beta * Q + lambda1 * V_N * (beta + Q)) / (beta * Q * V_N * (lambda1 - lambda2)) * exp(lambda2 * t)
    C_N <- term1 + term2 - term3
    return(C_N)
  }
  
  # 公式6-4: 時間濃度解 - 遠場
  calc_C_F_time <- function(t, G, beta, Q, V_N, V_F, lambda1, lambda2) {
    term1 <- G / Q
    term2 <- G * (lambda1 * V_N + beta) / beta * 
      (beta * Q + lambda2 * V_N * (beta + Q)) / (beta * Q * V_N * (lambda1 - lambda2)) * exp(lambda1 * t)
    term3 <- G * (lambda2 * V_N + beta) / beta * 
      (beta * Q + lambda1 * V_N * (beta + Q)) / (beta * Q * V_N * (lambda1 - lambda2)) * exp(lambda2 * t)
    C_F <- term1 + term2 - term3
    return(C_F)
  }
  
  # 公式6-5: 穩態濃度 - 近場
  calc_C_N_ss <- function(G, beta, Q) {
    C_N_SS <- G / Q + G / beta
    return(C_N_SS)
  }
  
  # 公式6-6: 穩態濃度 - 遠場
  calc_C_F_ss <- function(G, Q) {
    C_F_SS <- G / Q
    return(C_F_SS)
  }
  
  # 公式6-8: 濃度衰減 - 近場
  calc_C_N_decay <- function(t, C_N0, C_F0, beta, V_N, lambda1, lambda2) {
    term1 <- (beta * (C_F0 - C_N0) - lambda2 * V_N * C_N0) / (V_N * (lambda1 - lambda2)) * exp(lambda1 * t)
    term2 <- (beta * (C_N0 - C_F0) + lambda1 * V_N * C_N0) / (V_N * (lambda1 - lambda2)) * exp(lambda2 * t)
    C_N <- term1 + term2
    return(C_N)
  }
  
  # 公式6-9: 濃度衰減 - 遠場
  calc_C_F_decay <- function(t, C_N0, C_F0, beta, Q, V_N, V_F, lambda1, lambda2) {
    term1 <- ((lambda1 * V_N + beta) * (beta * (C_F0 - C_N0) - lambda2 * V_N * C_N0)) / 
      (beta * V_N * (lambda1 - lambda2)) * exp(lambda1 * t)
    term2 <- ((lambda2 * V_N + beta) * (beta * (C_N0 - C_F0) + lambda1 * V_N * C_N0)) / 
      (beta * V_N * (lambda1 - lambda2)) * exp(lambda2 * t)
    C_F <- term1 + term2
    return(C_F)
  }
  
  # 計算暴露風險
  calc_risk <- function(conc_TWA, conc_STEL, TWA_lim, STEL_lim, IDLH_lim, OEL_lim) {
    risk_short <- ifelse(conc_STEL > IDLH_lim, "等級4 (IDLH超標)",
                         ifelse(conc_STEL > STEL_lim, "等級4 (STEL超標)", "安全"))
    risk_long <- ifelse(conc_TWA < OEL_lim * 0.1, "等級1",
                        ifelse(conc_TWA < OEL_lim * 0.5, "等級2",
                               ifelse(conc_TWA < OEL_lim, "等級3", "等級4")))
    return(list(short = risk_short, long = risk_long))
  }
  
  # 主計算函數
  result_data <- eventReactive(input$calc, {
    G <- input$G
    Q <- input$Q
    beta <- input$beta
    V_N <- input$V_N
    V_F <- input$V_F
    t_work <- input$t_work
    t_decay_start <- input$t_decay_start
    t_total <- input$t_total
    
    lambdas <- calc_lambda(beta, Q, V_N, V_F)
    lambda1 <- lambdas$lambda1
    lambda2 <- lambdas$lambda2
    
    if (input$func == "時間濃度解") {
      # 完整模擬：排放階段 + 衰減階段
      t_seq_emission <- seq(0, t_decay_start, by = 1)
      C_N_emission <- sapply(t_seq_emission, function(t) 
        calc_C_N_time(t, G, beta, Q, V_N, V_F, lambda1, lambda2))
      C_F_emission <- sapply(t_seq_emission, function(t) 
        calc_C_F_time(t, G, beta, Q, V_N, V_F, lambda1, lambda2))
      
      # 衰減階段（從衰減開始時間起）
      t_seq_decay <- seq(t_decay_start + 1, t_total, by = 1)
      C_N0 <- tail(C_N_emission, 1)
      C_F0 <- tail(C_F_emission, 1)
      
      C_N_decay <- sapply(t_seq_decay - t_decay_start, function(t) 
        calc_C_N_decay(t, C_N0, C_F0, beta, V_N, lambda1, lambda2))
      C_F_decay <- sapply(t_seq_decay - t_decay_start, function(t) 
        calc_C_F_decay(t, C_N0, C_F0, beta, Q, V_N, V_F, lambda1, lambda2))
      
      # 合併數據
      df_emission <- data.frame(
        Time = t_seq_emission,
        C_N = C_N_emission,
        C_F = C_F_emission,
        Phase = "排放階段"
      )
      
      df_decay <- data.frame(
        Time = t_seq_decay,
        C_N = C_N_decay,
        C_F = C_F_decay,
        Phase = "衰減階段"
      )
      
      df <- rbind(df_emission, df_decay)
      
      # 計算暴露指標（僅考慮排放階段）
      TWA_N <- mean(C_N_emission) * min(t_decay_start, 480) / 480
      if(length(C_N_emission) > 14) {
        STEL_N <- max(sapply(1:(length(C_N_emission)-14), function(i) mean(C_N_emission[i:(i+14)])))
      } else {
        STEL_N <- mean(C_N_emission)
      }
      
    } else if (input$func == "穩態濃度") {
      # 使用公式6-5和6-6
      C_N_SS <- calc_C_N_ss(G, beta, Q)
      C_F_SS <- calc_C_F_ss(G, Q)
      
      df <- data.frame(
        Time = 0,
        C_N = C_N_SS,
        C_F = C_F_SS,
        Phase = "穩態"
      )
      
      TWA_N <- C_N_SS
      STEL_N <- C_N_SS
      
    } else if (input$func == "濃度衰減") {
      # 直接從衰減開始時間計算（假設初始濃度為穩態）
      C_N0 <- calc_C_N_ss(G, beta, Q)
      C_F0 <- calc_C_F_ss(G, Q)
      
      t_seq <- seq(t_decay_start, t_total, by = 1)
      C_N <- sapply(t_seq - t_decay_start, function(t) 
        calc_C_N_decay(t, C_N0, C_F0, beta, V_N, lambda1, lambda2))
      C_F <- sapply(t_seq - t_decay_start, function(t) 
        calc_C_F_decay(t, C_N0, C_F0, beta, Q, V_N, V_F, lambda1, lambda2))
      
      df <- data.frame(
        Time = t_seq,
        C_N = C_N,
        C_F = C_F,
        Phase = "衰減階段"
      )
      
      TWA_N <- mean(C_N) * (t_total - t_decay_start) / 480
      if(length(C_N) > 14) {
        STEL_N <- max(sapply(1:(length(C_N)-14), function(i) mean(C_N[i:(i+14)])))
      } else {
        STEL_N <- mean(C_N)
      }
    }
    
    risk <- calc_risk(TWA_N, STEL_N, input$TWA, input$STEL, input$IDLH, input$OEL)
    
    list(df = df, TWA_N = TWA_N, STEL_N = STEL_N, risk = risk, 
         t_decay_start = t_decay_start, t_total = t_total)
  })
  
  # 顯示計算結果
  output$result_text <- renderPrint({
    data <- result_data()
    cat("=== 計算結果 ===\n")
    cat("近場 TWA 濃度 (8hr):", round(data$TWA_N, 2), "mg/m³\n")
    cat("近場 STEL 濃度 (15min):", round(data$STEL_N, 2), "mg/m³\n")
    cat("\n=== 風險評估 ===\n")
    cat("短期暴露風險:", data$risk$short, "\n")
    cat("長期暴露風險:", data$risk$long, "\n")
    
    # 顯示使用的公式
    cat("\n=== 使用公式 ===\n")
    if (input$func == "時間濃度解") {
      cat("近場濃度: 公式6-3 (排放階段), 公式6-8 (衰減階段)\n")
      cat("遠場濃度: 公式6-4 (排放階段), 公式6-9 (衰減階段)\n")
      cat("衰減開始時間:", input$t_decay_start, "分鐘\n")
    } else if (input$func == "穩態濃度") {
      cat("近場濃度: 公式6-5\n")
      cat("遠場濃度: 公式6-6\n")
    } else if (input$func == "濃度衰減") {
      cat("近場濃度: 公式6-8\n")
      cat("遠場濃度: 公式6-9\n")
      cat("衰減開始時間:", input$t_decay_start, "分鐘\n")
    }
    cat("總模擬時間:", input$t_total, "分鐘\n")
  })
  
  # 顯示風險分級表
  output$risk_table <- renderTable({
    data <- result_data()
    df_risk <- data.frame(
      類型 = c("IDLH", "STEL", "TWA", "OEL"),
      限值 = c(input$IDLH, input$STEL, input$TWA, input$OEL),
      計算值 = c(round(data$STEL_N, 2), round(data$STEL_N, 2), round(data$TWA_N, 2), round(data$TWA_N, 2)),
      風險等級 = c(
        ifelse(data$STEL_N > input$IDLH, "等級4", ifelse(data$STEL_N > input$IDLH * 0.5, "等級3", "安全")),
        ifelse(data$STEL_N > input$STEL, "等級4", ifelse(data$STEL_N > input$STEL * 0.5, "等級3", "安全")),
        ifelse(data$TWA_N > input$TWA, "等級4", ifelse(data$TWA_N > input$TWA * 0.5, "等級3", "安全")),
        ifelse(data$TWA_N > input$OEL, "等級4", ifelse(data$TWA_N > input$OEL * 0.5, "等級3", "安全"))
      )
    )
    df_risk
  })
  
  # 繪製濃度趨勢圖
  output$conc_plot <- renderPlot({
    data <- result_data()
    df <- data$df
    
    if (nrow(df) > 1) {
      p <- ggplot(df, aes(x = Time)) +
        geom_line(aes(y = C_N, color = "近場"), linewidth = 1) +
        geom_line(aes(y = C_F, color = "遠場"), linewidth = 1) +
        geom_vline(xintercept = data$t_decay_start, linetype = "dashed", color = "red", linewidth = 1) +
        annotate("text", x = data$t_decay_start, y = max(df$C_N, na.rm = TRUE) * 0.95, 
                 label = paste("衰減開始:", data$t_decay_start, "分鐘"), 
                 hjust = -0.1, color = "red") +
        labs(title = paste0("近場與遠場濃度隨時間變化 (", input$func, ")"),
             x = "時間 (分鐘)", y = "濃度 (mg/m³)") +
        scale_color_manual(values = c("近場" = "red", "遠場" = "blue")) +
        theme_minimal() +
        theme(legend.position = "bottom")
    } else {
      p <- ggplot() +
        geom_point(aes(x = 0, y = df$C_N, color = "近場"), size = 5) +
        geom_point(aes(x = 0, y = df$C_F, color = "遠場"), size = 5) +
        labs(title = paste0("穩態濃度 (", input$func, ")"),
             x = "", y = "濃度 (mg/m³)") +
        scale_color_manual(values = c("近場" = "red", "遠場" = "blue")) +
        theme_minimal() +
        theme(legend.position = "bottom",
              axis.text.x = element_blank())
    }
    
    p
  })
}

# 執行 Shiny App
shinyApp(ui, server)