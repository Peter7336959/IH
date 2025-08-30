# --- 1. 載入所需套件 ---
library(shiny)
library(ggplot2)
library(dplyr)
library(DT)

# --- 2. 使用者介面 (UI) 定義 ---
ui <- fluidPage(
  withMathJax(), # 啟用MathJax來渲染LaTeX公式
  
  # 應用程式標題
  titlePanel("二暴露區模式(Two-Zone Model)暴露評估與風險管理工具"),
  
  # 插入模型示意圖
  # 注意：請在app.R檔案的同層目錄下建立一個名為 'www' 的資料夾，並將 '2zone.png' 圖片放入其中
  div(
    img(src = "2zone.png", alt = "暴露空間模式示意圖", style = "max-width: 30%;"),
    style = "text-align: center;"
  ),
  
  hr(),
  
  # 側邊欄佈局
  sidebarLayout(
    # --- 2.1 輸入面板 (Sidebar Panel) ---
    sidebarPanel(
      width = 4,
      h4("參數輸入"),
      tabsetPanel(
        # --- 標籤頁1: 情境與化學物質 ---
        tabPanel("情境與化學物質",
                 br(),
                 selectInput("calc_mode", "計算功能選擇:",
                             choices = c("時間濃度解" = "time_dependent",
                                         "穩態濃度" = "steady_state",
                                         "濃度衰減(污染源停止逸散)" = "decay")),
                 textInput("chem_name", "化學物質名稱:", "例如：甲苯"),
                 numericInput("mw", "分子量 (g/mol):", value = 92.14, min = 1),
                 numericInput("G", "污染源排放速率 (G, mg/min):", value = 100, min = 0),
                 numericInput("t_dur", "作業開始後時間 (t, 分鐘):", value = 120, min = 1),
                 conditionalPanel(
                   condition = "input.calc_mode == 'decay'",
                   numericInput("CN0", "近場初始濃度 (CN0, mg/m³):", value = 0, min = 0),
                   numericInput("CF0", "遠場初始濃度 (CF0, mg/m³):", value = 0, min = 0)
                 )
        ),
        # --- 標籤頁2: 物理環境 ---
        tabPanel("物理環境",
                 br(),
                 numericInput("V_total", "室內總體積 (m³):", value = 300, min = 1),
                 numericInput("r", "人員近場呼吸域半徑 (r, m):", value = 0.5, min = 0.1, step = 0.1),
                 numericInput("S", "近場空氣流速 (S, m/min):", value = 10, min = 0.1),
                 numericInput("Q", "房間通風換氣速率 (Q, m³/min):", value = 30, min = 0.1)
        ),
        # --- 標籤頁3: 暴露限值 ---
        tabPanel("暴露限值設定",
                 br(),
                 p("請輸入暴露限值，若無資料可留白或填0。"),
                 fluidRow(
                   column(6, numericInput("twa_limit", "TWA (8小時):", value = 50)),
                   column(6, radioButtons("twa_unit", "單位:", choices = c("ppm", "mg/m³"), selected = "ppm", inline = TRUE))
                 ),
                 fluidRow(
                   column(6, numericInput("stel_limit", "STEL (15分鐘):", value = 100)),
                   column(6, radioButtons("stel_unit", "單位:", choices = c("ppm", "mg/m³"), selected = "ppm", inline = TRUE))
                 ),
                 fluidRow(
                   column(6, numericInput("idlh_limit", "IDLH:", value = 500)),
                   column(6, radioButtons("idlh_unit", "單位:", choices = c("ppm", "mg/m³"), selected = "ppm", inline = TRUE))
                 ),
                 fluidRow(
                   column(6, numericInput("oel_limit", "其他OEL:", value = 50)),
                   column(6, radioButtons("oel_unit", "單位:", choices = c("ppm", "mg/m³"), selected = "ppm", inline = TRUE))
                 )
        )
      ),
      hr(),
      actionButton("calculate", "開始計算", class = "btn-primary", style = "width: 100%;"),
      br(), br(),
      # 簡單的列印按鈕
      tags$button("列印本頁", onclick = "window.print();", class = "btn btn-info", style = "width: 100%;")
    ),
    
    # --- 2.2 輸出面板 (Main Panel) ---
    mainPanel(
      width = 8,
      tabsetPanel(
        # --- 輸出標籤頁A: 風險評估儀表板 ---
        tabPanel("風險評估儀表板",
                 h4("衍生參數計算結果"),
                 verbatimTextOutput("derived_params_output"),
                 h4("風險評估摘要"),
                 DTOutput("risk_table"),
                 h4("綜合風險評級與控制建議"),
                 uiOutput("risk_summary_recommendation")
        ),
        # --- 輸出標籤頁B: 濃度時間趨勢圖 ---
        tabPanel("濃度時間趨勢圖",
                 plotOutput("concentration_plot", height = "500px")
        ),
        # --- 輸出標籤頁C: 模型公式參考 ---
        tabPanel("模型公式參考",
                 uiOutput("formula_display")
        )
      )
    )
  )
)

# --- 3. 伺服器端 (Server) 邏輯定義 ---
server <- function(input, output) {
  
  # --- 3.1 輔助函數 ---
  # 單位轉換函數: ppm -> mg/m3
  ppm_to_mgm3 <- function(ppm, mw) {
    # 修正：使用 |
    
     避免語法錯誤
    if (is.na(ppm) |
        
         is.na(mw) |
         mw <= 0) return(NA)
    return((ppm * mw) / 24.45)
  }
  
  # --- 3.2 核心計算邏輯 (使用eventReactive) ---
  results <- eventReactive(input$calculate, {
    
    # --- 3.2.1 收集與計算衍生參數 ---
    r <- input$r
    S <- input$S
    V_total <- input$V_total
    Q <- input$Q
    G <- input$G
    
    FSA <- 2 * pi * r^2
    VN <- (2/3) * pi * r^3
    VF <- V_total - VN
    beta <- 0.5 * FSA * S
    
    # 檢查VF是否為正
    if (VF <= 0) {
      showNotification("錯誤：遠場體積(VF)必須大於0。請檢查室內總體積與近場半徑。", type = "error")
      return(NULL)
    }
    
    # 計算 lambda 1 & 2
    term1 <- -((beta * VF + VN * (beta + Q)) / (VN * VF))
    term2_sqrt_base <- ((beta * VF + VN * (beta + Q)) / (VN * VF))^2 - 4 * ((beta * Q) / (VN * VF))
    
    if (term2_sqrt_base < 0) {
      showNotification("錯誤：Lambda計算出現問題(根號內為負)，請檢查輸入參數。", type = "error")
      return(NULL)
    }
    
    lambda1 <- 0.5 * (term1 + sqrt(term2_sqrt_base))
    lambda2 <- 0.5 * (term1 - sqrt(term2_sqrt_base))
    
    # 修正：增加數值穩定性檢查，避免 lambda1 和 lambda2 過於接近導致除以零
    if (abs(lambda1 - lambda2) < 1e-9) {
      showNotification("錯誤：模型參數導致數學不穩定(λ₁ ≈ λ₂)，請調整輸入值。", type = "error")
      return(NULL)
    }
    
    # --- 3.2.2 濃度計算 ---
    time_seq <- 1:input$t_dur
    conc_data <- data.frame(Time = time_seq, CN = 0, CF = 0)
    
    if (input$calc_mode == "time_dependent") {
      # 時間濃度解 (公式 7-1, 7-2)
      C_N_ss = G/Q + G/beta
      
      term_N1 = G / (lambda1 - lambda2) * ( (beta*Q + lambda2*VN*(beta+Q)) / (beta*Q*VN) )
      term_N2 = G / (lambda1 - lambda2) * ( (beta*Q + lambda1*VN*(beta+Q)) / (beta*Q*VN) )
      
      term_F1 = G / (lambda1 - lambda2) * ( (lambda1*VN+beta)/beta ) * ( (beta*Q + lambda2*VN*(beta+Q)) / (beta*Q*VN) )
      term_F2 = G / (lambda1 - lambda2) * ( (lambda2*VN+beta)/beta ) * ( (beta*Q + lambda1*VN*(beta+Q)) / (beta*Q*VN) )
      
      conc_data$CN <- C_N_ss + term_N1 * exp(lambda1 * time_seq) - term_N2 * exp(lambda2 * time_seq)
      conc_data$CF <- G/Q + term_F1 * exp(lambda1 * time_seq) - term_F2 * exp(lambda2 * time_seq)
      
    } else if (input$calc_mode == "steady_state") {
      # 穩態濃度 (公式 7-3, 7-4)
      conc_data$CN <- G / Q + G / beta
      conc_data$CF <- G / Q
      
    } else if (input$calc_mode == "decay") {
      # 濃度衰減 (公式 7-5, 7-6)
      CN0 <- input$CN0
      CF0 <- input$CF0
      
      term_N1 = (beta * (CF0 - CN0) - lambda2 * VN * CN0) / (VN * (lambda1 - lambda2))
      term_N2 = (beta * (CN0 - CF0) + lambda1 * VN * CN0) / (VN * (lambda1 - lambda2))
      
      term_F1 = ((lambda1*VN+beta)/beta) * term_N1
      term_F2 = ((lambda2*VN+beta)/beta) * term_N2
      
      conc_data$CN <- term_N1 * exp(lambda1 * time_seq) + term_N2 * exp(lambda2 * time_seq)
      conc_data$CF <- term_F1 * exp(lambda1 * time_seq) + term_F2 * exp(lambda2 * time_seq)
    }
    
    # --- 3.2.3 暴露指標計算 ---
    twa_calc <- sum(conc_data$CN) / 480
    
    # 計算15分鐘移動平均的最大值
    if (nrow(conc_data) >= 15) {
      stel_calc <- max(sapply(1:(nrow(conc_data) - 14), function(i) {
        mean(conc_data$CN[i:(i + 14)])
      }))
    } else {
      stel_calc <- mean(conc_data$CN) # 如果時間不足15分鐘，則取總平均
    }
    
    idlh_metric <- max(conc_data$CN, na.rm = TRUE)
    
    # --- 3.2.4 單位標準化與風險評估 ---
    mw <- input$mw
    limits <- list(
      TWA = list(val = input$twa_limit, unit = input$twa_unit),
      STEL = list(val = input$stel_limit, unit = input$stel_unit),
      IDLH = list(val = input$idlh_limit, unit = input$idlh_unit),
      OEL = list(val = input$oel_limit, unit = input$oel_unit)
    )
    
    risk_df <- data.frame()
    
    for (name in names(limits)) {
      limit_val <- limits[[name]]$val
      limit_unit <- limits[[name]]$unit
      
      if (!is.na(limit_val) && limit_val > 0) {
        limit_mgm3 <- ifelse(limit_unit == "ppm", ppm_to_mgm3(limit_val, mw), limit_val)
        
        assessed_conc <- switch(name,
                                TWA = twa_calc,
                                OEL = twa_calc,
                                STEL = stel_calc,
                                IDLH = idlh_metric)
        
        ratio <- assessed_conc / limit_mgm3
        
        # 風險分級邏輯
        level <- 1
        status <- "可接受"
        if (name %in% c("IDLH", "STEL")) {
          if (ratio > 1) {
            level <- 4
            status <- "超過限值"
          }
        } else { # TWA, OEL
          if (ratio > 1) {
            level <- 4
            status <- "超過限值"
          } else if (ratio > 0.5) {
            level <- 3
            status <- "警戒水平"
          } else if (ratio > 0.1) {
            level <- 2
            status <- "需管理"
          }
        }
        
        risk_df <- rbind(risk_df, data.frame(
          "暴露限值項目" = name,
          "限值 (輸入單位)" = paste(limit_val, limit_unit),
          "限值 (mg/m³)" = round(limit_mgm3, 2),
          "評估暴露濃度 (mg/m³)" = round(assessed_conc, 2),
          "風險比 (評估/限值)" = round(ratio, 2),
          "風險分級 (1-4)" = level,
          "狀態" = status,
          check.names = FALSE
        ))
      }
    }
    
    # --- 3.2.5 返回所有結果 ---
    # 建立一個輔助函數，安全地從風險表中提取單一純量值給繪圖函數
    get_limit_mgm3 <- function(df, name) {
      val <- df$`限值 (mg/m³)`[df$`暴露限值項目` == name]
      if (length(val) == 1) return(val) else return(NA)
    }
    
    return(list(
      derived_params = list(
        FSA = FSA, VN = VN, VF = VF, beta = beta, lambda1 = lambda1, lambda2 = lambda2,
        CNF = tail(conc_data$CN, 1), # 新增：最終近場濃度
        CFF = tail(conc_data$CF, 1)  # 新增：最終遠場濃度
      ),
      conc_data = conc_data,
      risk_df = risk_df,
      limits_mgm3 = list(
        TWA = get_limit_mgm3(risk_df, "TWA"),
        STEL = get_limit_mgm3(risk_df, "STEL"),
        IDLH = get_limit_mgm3(risk_df, "IDLH")
      )
    ))
  })
  
  # --- 3.3 輸出渲染 ---
  
  # 渲染衍生參數
  output$derived_params_output <- renderPrint({
    res <- results()
    if (is.null(res)) return("請點擊'開始計算'以生成結果。")
    
    cat(paste("近場自由表面積 (FSA):", round(res$derived_params$FSA, 2), "m²\n"))
    cat(paste("近場體積 (VN):", round(res$derived_params$VN, 2), "m³\n"))
    cat(paste("遠場體積 (VF):", round(res$derived_params$VF, 2), "m³\n"))
    cat(paste("區域間空氣交換速率 (β):", round(res$derived_params$beta, 2), "m³/min\n"))
    cat(paste("特徵值 (λ₁):", round(res$derived_params$lambda1, 4), "min⁻¹\n"))
    cat(paste("特徵值 (λ₂):", round(res$derived_params$lambda2, 4), "min⁻¹\n"))
    cat(paste("最終近場濃度 (CNF):", round(res$derived_params$CNF, 2), "mg/m³\n"))
    cat(paste("最終遠場濃度 (CFF):", round(res$derived_params$CFF, 2), "mg/m³\n"))
  })
  
  # 渲染風險評估表
  output$risk_table <- renderDT({
    res <- results()
    # 修正：使用 |
    
     避免語法錯誤
    if (is.null(res) |
        
         nrow(res$risk_df) == 0) {
      return(datatable(data.frame(訊息 = "無風險評估結果可顯示。")))
    }
    datatable(res$risk_df, options = list(pageLength = 5, dom = 't'), rownames = FALSE)
  })
  
  # 渲染綜合風險建議
  output$risk_summary_recommendation <- renderUI({
    res <- results()
    # 修正：使用 |
    
     避免語法錯誤
    if (is.null(res) |
        
         nrow(res$risk_df) == 0) return(p("無結果可供評級。"))
    
    max_level <- max(res$risk_df$`風險分級 (1-4)`)
    
    level_text <- paste0("<h4>綜合風險等級: ", max_level, "</h4>")
    
    recommendation_text <- switch(as.character(max_level),
                                  "4" = "<b>高度風險</b>: 立即停止作業。暴露濃度已超過短期或八小時暴露限值，可能對生命健康造成立即或嚴重威脅。<br><b>建議措施</b>: 必須採取工程控制或源頭替代/消除措施。例如：安裝局部排氣通風系統、更換為低危害性化學品、或改變製程完全消除逸散源。",
                                  "3" = "<b>中度風險</b>: 應立即採取行動降低暴露。暴露濃度已超過行動水平(50% OEL)。<br><b>建議措施</b>: 優先考慮行政管理控制並檢討工程控制有效性。例如：減少工作時間、增加休息頻率、改善作業程序、加強人員訓練。",
                                  "2" = "<b>低度風險</b>: 需進行風險管理。暴露濃度低於行動水平，但仍需注意。<br><b>建議措施</b>: 確保現有控制措施有效。例如：確認個人防護具(PPE)的正確使用與維護、定期環境監測、員工健康管理。",
                                  "1" = "<b>可接受風險</b>: 維持現行做法。暴露濃度在可接受範圍內。<br><b>建議措施</b>: 持續維持良好的作業習慣與現有控制措施。"
    )
    
    HTML(paste(level_text, p(HTML(recommendation_text))))
  })
  
  # 渲染濃度趨勢圖
  output$concentration_plot <- renderPlot({
    res <- results()
    if (is.null(res)) return(NULL)
    
    p <- ggplot(res$conc_data, aes(x = Time)) +
      geom_line(aes(y = CN, color = "近場濃度 (CN)"), size = 1.2) +
      geom_line(aes(y = CF, color = "遠場濃度 (CF)"), size = 1.2) +
      labs(title = "近場與遠場濃度隨時間變化趨勢圖",
           x = "時間 (分鐘)",
           y = "濃度 (mg/m³)",
           color = "圖例") +
      theme_minimal(base_size = 16) +
      scale_color_manual(values = c("近場濃度 (CN)" = "red", "遠場濃度 (CF)" = "blue"))
    
    # 添加暴露限值水平線
    if(!is.na(res$limits_mgm3$TWA)) {
      p <- p + geom_hline(yintercept = res$limits_mgm3$TWA, linetype = "dashed", color = "darkgreen") +
        annotate("text", x = max(res$conc_data$Time)*0.9, y = res$limits_mgm3$TWA, label = "TWA Limit", vjust = -0.5, color = "darkgreen")
    }
    if(!is.na(res$limits_mgm3$STEL)) {
      p <- p + geom_hline(yintercept = res$limits_mgm3$STEL, linetype = "dashed", color = "orange") +
        annotate("text", x = max(res$conc_data$Time)*0.9, y = res$limits_mgm3$STEL, label = "STEL Limit", vjust = -0.5, color = "orange")
    }
    if(!is.na(res$limits_mgm3$IDLH)) {
      p <- p + geom_hline(yintercept = res$limits_mgm3$IDLH, linetype = "dashed", color = "purple") +
        annotate("text", x = max(res$conc_data$Time)*0.9, y = res$limits_mgm3$IDLH, label = "IDLH Limit", vjust = -0.5, color = "purple")
    }
    
    print(p)
  })
  
  # 渲染公式
  output$formula_display <- renderUI({
    mode <- input$calc_mode
    
    formula_text <- switch(mode,
                           "time_dependent" = "<h4>時間濃度解 (公式 7-1 & 7-2)</h4>
        <p><b>近場濃度 \\(C_N(t)\\):</b></p>
        $$C_{N}(t) = \\frac{G}{Q} + \\frac{G}{\\beta} + G[\\frac{\\beta Q+\\lambda_{2}V_{N}(\\beta+Q)}{\\beta QV_{N}(\\lambda_{1}-\\lambda_{2})}]e^{\\lambda_{1}t} - G[\\frac{\\beta Q+\\lambda_{1}V_{N}(\\beta+Q)}{\\beta QV_{N}(\\lambda_{1}-\\lambda_{2})}]e^{\\lambda_{2}t}$$
        <p><b>遠場濃度 \\(C_F(t)\\):</b></p>
        $$C_{F}(t) = \\frac{G}{Q} + G[\\frac{\\lambda_{1}V_{N}+\\beta}{\\beta}][\\frac{\\beta Q+\\lambda_{2}V_{N}(\\beta+Q)}{\\beta QV_{N}(\\lambda_{1}-\\lambda_{2})}]e^{\\lambda_{1}t} - G[\\frac{\\lambda_{2}V_{N}+\\beta}{\\beta}][\\frac{\\beta Q+\\lambda_{1}V_{N}(\\beta+Q)}{\\beta QV_{N}(\\lambda_{1}-\\lambda_{2})}]e^{\\lambda_{2}t}$$
        ",
                           "steady_state" = "<h4>穩態濃度 (公式 7-3 & 7-4)</h4>
        <p><b>近場穩態濃度 \\(C_{N,SS}\\):</b></p>
        $$C_{N,SS} = \\frac{G}{Q} + \\frac{G}{\\beta}$$
        <p><b>遠場穩態濃度 \\(C_{F,SS}\\):</b></p>
        $$C_{F,SS} = \\frac{G}{Q}$$
        ",
                           "decay" = "<h4>濃度衰減 (公式 7-5 & 7-6)</h4>
        <p><b>近場濃度衰減 \\(C_N(t)\\):</b></p>
        $$C_{N}(t) = \\frac{\\beta(C_{F_{0}}-C_{N_{0}})-\\lambda_{2}V_{N}C_{N_{0}}}{V_{N}(\\lambda_{1}-\\lambda_{2})}e^{\\lambda_{1}t} + \\frac{\\beta(C_{N_{0}}-C_{F_{0}})+\\lambda_{1}V_{N}C_{N_{0}}}{V_{N}(\\lambda_{1}-\\lambda_{2})}e^{\\lambda_{2}t}$$
        <p><b>遠場濃度衰減 \\(C_F(t)\\):</b></p>
        $$C_{F}(t) = [\\frac{\\lambda_{1}V_{N}+\\beta}{\\beta}] \\frac{\\beta(C_{F_{0}}-C_{N_{0}})-\\lambda_{2}V_{N}C_{N_{0}}}{V_{N}(\\lambda_{1}-\\lambda_{2})}e^{\\lambda_{1}t} + [\\frac{\\lambda_{2}V_{N}+\\beta}{\\beta}] \\frac{\\beta(C_{N_{0}}-C_{F_{0}})+\\lambda_{1}V_{N}C_{N_{0}}}{V_{N}(\\lambda_{1}-\\lambda_{2})}e^{\\lambda_{2}t}$$
        "
    )
    withMathJax(HTML(formula_text))
  })
}

# --- 4. 執行應用程式 ---
shinyApp(ui = ui, server = server)