"""
HTML templates and CSS styling for benchmark report generation.

This module contains all HTML/CSS templates used to build the interactive
HTML report with embedded visualizations.
"""

from typing import Dict, Any


def get_html_template() -> str:
    """Return the HTML template with embedded CSS and inline Plotly JS.

    The Plotly JS bundle is embedded directly so the resulting HTML works
    without an internet connection.

    Returns:
        HTML template string with {placeholders} for content injection.
    """
    from plotly.offline import get_plotlyjs as _get_plotlyjs
    plotlyjs_tag = f'<script type="text/javascript">{_get_plotlyjs()}</script>'

    template = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{title}</title>
    __PLOTLYJS__
    <style>
        /* Base Styles */
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}
        
        body {{
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Arial, sans-serif;
            line-height: 1.6;
            color: #333;
            background: #f5f7fa;
        }}
        
        /* Header */
        .header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 2rem;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
        }}
        
        .header h1 {{
            font-size: 2.5rem;
            margin-bottom: 0.5rem;
        }}
        
        /* Navigation */
        .nav {{
            background: white;
            padding: 1rem 2rem;
            box-shadow: 0 2px 4px rgba(0,0,0,0.05);
            position: sticky;
            top: 0;
            z-index: 100;
        }}
        
        .nav a {{
            color: #667eea;
            text-decoration: none;
            margin-right: 2rem;
            font-weight: 500;
            transition: all 0.3s;
        }}
        
        .nav a:hover {{
            color: #764ba2;
        }}
        
        /* Main Container */
        .container {{
            max-width: 100%;
            margin: 2rem auto;
            padding: 0 2rem;
        }}
        
        /* Sections */
        .section {{
            background: white;
            padding: 2rem;
            margin-bottom: 2rem;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.05);
        }}
        
        .section h2 {{
            color: #667eea;
            margin-bottom: 1.5rem;
            padding-bottom: 0.5rem;
            border-bottom: 2px solid #667eea;
        }}
        
        /* Tables */
        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 1rem 0;
        }}
        
        th {{
            background: #667eea;
            color: white;
            padding: 0.75rem;
            text-align: left;
            font-weight: 600;
        }}
        
        td {{
            padding: 0.75rem;
            border-bottom: 1px solid #e0e0e0;
        }}
        
        tr:hover {{
            background: #f9f9f9;
        }}
        
        /* Tabs */
        .tabs {{
            display: flex;
            border-bottom: 2px solid #e0e0e0;
            margin-bottom: 1.5rem;
        }}
        
        .tab {{
            padding: 0.75rem 1.5rem;
            cursor: pointer;
            border: none;
            background: none;
            font-size: 1rem;
            color: #666;
            transition: all 0.3s;
        }}
        
        .tab:hover {{
            background: #f0f0f0;
        }}
        
        .tab.active {{
            color: #667eea;
            border-bottom: 3px solid #667eea;
            font-weight: 600;
        }}
        
        .tab-content {{
            display: none;
        }}
        
        .tab-content.active {{
            display: block;
        }}
        
        /* Plots */
        .plot-container {{
            margin: 2rem 0;
            width: 100%;
            overflow-x: visible;
        }}
        
        /* Force Plotly divs to fill container */
        .plot-container > div,
        .plot-box > div,
        .js-plotly-plot {{
            width: 100% !important;
        }}

        /* Flex layout for plot + controls */
        .plot-flex {{
            display: flex;
            gap: 1rem;
            align-items: flex-start;
            width: 100%;
        }}

        .plot-box {{
            flex: 1 1 auto;
            min-width: 0; /* allow flex child to shrink so plotly's div can be responsive */
            width: 100%;
        }}

        .plot-controls {{
            flex: 0 0 200px;
            width: 200px;
            max-height: 600px;
            overflow-y: auto;
            overflow-x: hidden;
            background: #fbfbfb;
            border: 1px solid #ebebeb;
            padding: 0.75rem;
            border-radius: 6px;
            box-shadow: 0 1px 2px rgba(0,0,0,0.04);
            font-size: 13px;
        }}

        .plot-controls .controls-header {{ 
            margin-bottom: 0.5rem; 
            position: sticky;
            top: 0;
            background: #fbfbfb;
            padding-bottom: 0.25rem;
            border-bottom: 1px solid #e0e0e0;
        }}
        
        .plot-controls .controls-section {{ 
            margin-bottom: 0.75rem; 
        }}
        
        .plot-controls label {{
            display: block;
            margin: 0.25rem 0;
            cursor: pointer;
        }}
        
        .plot-controls label:hover {{
            background: #f0f0f0;
            padding-left: 2px;
        }}
        
        /* Info Box */
        .info-box {{
            background: #e3f2fd;
            border-left: 4px solid #2196f3;
            padding: 1rem;
            margin: 1rem 0;
            border-radius: 4px;
        }}
        
        /* Footer */
        .footer {{
            text-align: center;
            padding: 2rem;
            color: #666;
            font-size: 0.9rem;
        }}
    </style>
    <script>
        function switchTab(tabName, event) {{
            // Hide all tab contents
            document.querySelectorAll('.tab-content').forEach(content => {{
                content.classList.remove('active');
            }});
            
            // Remove active from all tabs
            document.querySelectorAll('.tab').forEach(tab => {{
                tab.classList.remove('active');
            }});
            
            // Show selected tab content
            document.getElementById(tabName).classList.add('active');
            if (event && event.target) {{
                event.target.classList.add('active');
            }}
            
            // Force resize all Plotly plots in the newly visible tab
            setTimeout(function() {{
                var activeTab = document.getElementById(tabName);
                if (activeTab) {{
                    var plots = activeTab.querySelectorAll('.js-plotly-plot');
                    plots.forEach(function(plot) {{
                        Plotly.Plots.resize(plot);
                    }});
                }}
            }}, 50);
        }}
        
        // Auto-resize all plots when page loads and when window resizes
        window.addEventListener('load', function() {{
            setTimeout(function() {{
                var allPlots = document.querySelectorAll('.js-plotly-plot');
                allPlots.forEach(function(plot) {{
                    Plotly.Plots.resize(plot);
                }});
            }}, 200);
        }});
        
        window.addEventListener('resize', function() {{
            var allPlots = document.querySelectorAll('.js-plotly-plot');
            allPlots.forEach(function(plot) {{
                Plotly.Plots.resize(plot);
            }});
        }});
    </script>
</head>
<body>
    <div class="header">
        <h1>{title}</h1>
        <p>{subtitle}</p>
    </div>
    
    <div class="nav">
        <a href="#summary">Summary</a>
        <a href="#experiment">Experiment Details</a>
        <a href="#visualizations">Visualizations</a>
    </div>
    
    <div class="container">
        {content}
    </div>
    
    <div class="footer">
        Generated by SigFitTest Workflow | {timestamp}
    </div>
</body>
</html>
"""
    return template.replace('__PLOTLYJS__', plotlyjs_tag)


def format_tool_name(tool: str) -> str:
    """Format tool name for display.
    
    Args:
        tool: Tool name (e.g., 'musical', 'spa')
        
    Returns:
        Formatted tool name (e.g., 'MuSiCal', 'SPA')
    """
    if tool.lower() == 'musical':
        return 'MuSiCal'
    elif tool.lower() == 'spa':
        return tool.upper()
    else:
        return tool


def create_info_box(title: str, content: str, color: str = '#2196f3') -> str:
    """Create an informational box with custom styling.
    
    Args:
        title: Box title
        content: HTML content to display
        color: Border color (default: blue)
        
    Returns:
        HTML string for info box
    """
    return f"""
    <div class="info-box" style="border-left-color: {color};">
        <h4 style="margin-bottom: 0.5rem;">{title}</h4>
        {content}
    </div>
    """


def create_feature_metric_scatter_widget(
        features_tsv_path: str,
        per_sample_metrics_tsv: str,
        widget_id: str = 'feature_metric_widget',
        data_df = None
) -> str:
        """Return HTML + JS for interactive boxplot: feature binned low/high vs metric distribution.

        Args:
                features_tsv_path: path to TSV file with scalar features per sample (DEPRECATED if data_df provided)
                per_sample_metrics_tsv: path to TSV with per-sample metrics (DEPRECATED if data_df provided)
                widget_id: unique DOM id prefix for this widget
                data_df: pandas DataFrame with all data (preferred over file paths to avoid CORS issues)

        The function produces a small interactive UI where the user selects:
            - feature (to bin into low/high)
            - metric (Y-axis for boxplot)
            - mutation burden filter

        Shows separate boxplots for each tool, split by low vs high feature values.
        Output: HTML string to embed in report.
        """
        import json

        html = f"""
<div style="background: white; padding: 1.5rem; border-radius: 8px; box-shadow: 0 2px 8px rgba(0,0,0,0.1); margin: 1rem 0;">
    <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(250px, 1fr)); gap: 1rem; margin-bottom: 1.5rem;">
        <div>
            <label style="display: block; font-weight: 600; color: #374151; margin-bottom: 0.5rem;">
                📊 Feature (binned Low/High)
            </label>
            <select id="{widget_id}_feature" style="width: 100%; padding: 0.5rem; border: 2px solid #e5e7eb; border-radius: 6px; font-size: 14px; background: white; cursor: pointer;">
                <option value="">Select feature...</option>
            </select>
        </div>
        <div>
            <label style="display: block; font-weight: 600; color: #374151; margin-bottom: 0.5rem;">
                📈 Metric (distribution)
            </label>
            <select id="{widget_id}_metric" style="width: 100%; padding: 0.5rem; border: 2px solid #e5e7eb; border-radius: 6px; font-size: 14px; background: white; cursor: pointer;">
                <option value="">Select metric...</option>
            </select>
        </div>
        <div>
            <label style="display: block; font-weight: 600; color: #374151; margin-bottom: 0.5rem;">
                🧬 Mutation Burden
            </label>
            <select id="{widget_id}_muts" style="width: 100%; padding: 0.5rem; border: 2px solid #e5e7eb; border-radius: 6px; font-size: 14px; background: white; cursor: pointer;">
                <option value="">All</option>
            </select>
        </div>
    </div>
    
    <div style="text-align: center; margin-bottom: 1rem;">
        <button id="{widget_id}_update" style="
            background: linear-gradient(135deg, #7c3aed 0%, #5b21b6 100%);
            color: white;
            border: none;
            padding: 0.75rem 2rem;
            border-radius: 8px;
            font-size: 16px;
            font-weight: 600;
            cursor: pointer;
            box-shadow: 0 4px 6px rgba(124, 58, 237, 0.3);
            transition: all 0.3s ease;
        " onmouseover="this.style.transform='translateY(-2px)'; this.style.boxShadow='0 6px 12px rgba(124, 58, 237, 0.4)';" 
           onmouseout="this.style.transform='translateY(0)'; this.style.boxShadow='0 4px 6px rgba(124, 58, 237, 0.3)';">
            🎯 Update Plot
        </button>
    </div>
    
    <div id="{widget_id}_plot" style="width: 100%; height: 600px; border: 2px solid #f3f4f6; border-radius: 8px; background: #fafafa;"></div>
    
    <div id="{widget_id}_stats" style="margin-top: 1rem; padding: 1rem; background: #f9fafb; border-radius: 6px; display: none;">
        <div style="display: flex; justify-content: center; align-items: center; gap: 2rem; font-size: 14px;">
            <div>
                <span style="color: #6b7280;">Median threshold:</span>
                <strong id="{widget_id}_threshold" style="color: #111827; margin-left: 0.5rem;">-</strong>
            </div>
            <div>
                <span style="color: #6b7280;">N (Low):</span>
                <strong id="{widget_id}_n_low" style="color: #111827; margin-left: 0.5rem;">-</strong>
            </div>
            <div>
                <span style="color: #6b7280;">N (High):</span>
                <strong id="{widget_id}_n_high" style="color: #111827; margin-left: 0.5rem;">-</strong>
            </div>
        </div>
    </div>
</div>

<script>
(function() {{
    const widgetId = '{widget_id}';
    
    // Feature columns (complexity features)
    const featureCols = ['highest_ENP_known', 'highest_ENP_unknown'];
    
    // Metric columns (performance metrics)
    const metricCols = ['TAE', 'F1', 'MCC', 'P', 'R', 'theoretical_min_TAE'];
    
    // Metadata columns
    const metadataCols = ['data_type', 'cancer_type', 'num_muts_condition', 'relevance', 
                          'tool', 'realization', 'sample', 'num_muts'];

    // Embedded data (avoids CORS issues with file:// protocol)
    const embeddedData = __DATA__;
    
    const features = embeddedData;
    const metrics = embeddedData;  // Same data source

    // Build feature list (only feature columns)
    const featureSelect = document.getElementById(widgetId + '_feature');
    console.log('Feature dropdown element:', featureSelect, 'Looking for ID:', widgetId + '_feature');
    console.log('Feature columns in data:', features.header.filter(h => featureCols.includes(h)));
    
    features.header.forEach(h => {{
        if (featureCols.includes(h)) {{
            const opt = document.createElement('option'); 
            opt.value = h; 
            opt.text = h; 
            featureSelect.appendChild(opt);
        }}
    }});

    // Build metric list (only metric columns)
    const metricSelect = document.getElementById(widgetId + '_metric');
    console.log('Metric dropdown element:', metricSelect, 'Looking for ID:', widgetId + '_metric');
    console.log('Metric columns in data:', metrics.header.filter(h => metricCols.includes(h)));
    
    metrics.header.forEach(h => {{
        if (metricCols.includes(h)) {{
            const opt = document.createElement('option'); 
            opt.value = h; 
                opt.text = h; 
                metricSelect.appendChild(opt);
            }}
        }});

        // Build mutation burdens list
        const toolIdx = metrics.header.indexOf('tool');
        const mutsIdx = metrics.header.indexOf('num_muts');
        const muts = new Set();
        if (mutsIdx >= 0) {{
            metrics.rows.forEach(r => {{
                if (r[mutsIdx]) muts.add(r[mutsIdx]);
            }});
        }}
        const mutsSelect = document.getElementById(widgetId + '_muts');
        Array.from(muts).sort((a,b)=>+a-+b).forEach(m => {{ 
            const opt = document.createElement('option'); 
            opt.value = m; 
            opt.text = m; 
            mutsSelect.appendChild(opt);
        }});

        // Build data structure for plotting
        const sampleIdx = metrics.header.indexOf('sample');
        const dataByToolMuts = new Map(); // tool_muts -> array of {{feature_vals, metric_vals, sample}}
        
        metrics.rows.forEach((row, idx) => {{
            const sample = row[sampleIdx];
            const tool = toolIdx >= 0 ? row[toolIdx] : '';
            const numMuts = mutsIdx >= 0 ? row[mutsIdx] : '';
            const key = tool + '_' + numMuts;
            
            // Get feature values for this row
            const featureVals = {{}};
            features.header.forEach((h, i) => {{
                if (featureCols.includes(h)) {{
                    featureVals[h] = row[i];  // Same row if same file
                }}
            }});
            
            // Get metric values for this row
            const metricVals = {{}};
            metrics.header.forEach((h, i) => {{
                if (metricCols.includes(h)) {{
                    metricVals[h] = row[i];
                }}
            }});
            
            if (!dataByToolMuts.has(key)) {{
                dataByToolMuts.set(key, []);
            }}
            dataByToolMuts.get(key).push({{ 
                sample, 
                tool, 
                numMuts, 
                features: featureVals, 
                metrics: metricVals 
            }});
        }});

        function updatePlot() {{
            const feat = featureSelect.value;
            const metric = metricSelect.value;
            
            if (!feat || !metric) {{
                document.getElementById('{widget_id}_plot').innerHTML = 
                    '<div style="display:flex; align-items:center; justify-content:center; height:100%; color:#9ca3af; font-size:18px;">📊 Select a feature and metric to visualize</div>';
                return;
            }}
            
            const selectedMuts = mutsSelect.value;

            // Collect all data points with feature and metric values
            const allData = [];
            
            dataByToolMuts.forEach((dataPoints, key) => {{
                dataPoints.forEach(point => {{
                    const numMuts = point.numMuts;
                    
                    // Filter by selected mutation burden
                    if (selectedMuts && numMuts !== selectedMuts) return;
                    
                    // Get feature and metric values
                    const featureVal = parseFloat(point.features[feat]);
                    const metricVal = parseFloat(point.metrics[metric]);
                    
                    // Skip if invalid
                    if (isNaN(featureVal) || isNaN(metricVal)) return;
                    
                    allData.push({{
                        tool: point.tool,
                        feature: featureVal,
                        metric: metricVal,
                        sample: point.sample,
                        numMuts: numMuts
                    }});
                }});
            }});
            
            if (allData.length === 0) {{
                document.getElementById('{widget_id}_plot').innerHTML = 
                    '<div style="display:flex; align-items:center; justify-content:center; height:100%; color:#ef4444; font-size:16px;">❌ No data available for selected options</div>';
                return;
            }}
            
            // Calculate median threshold for binning (proper median calculation)
            const featureValues = allData.map(d => d.feature);
            featureValues.sort((a,b) => a - b);
            let medianThreshold;
            if (featureValues.length % 2 === 0) {{
                // Even number: average of two middle values
                medianThreshold = (featureValues[featureValues.length/2 - 1] + featureValues[featureValues.length/2]) / 2;
            }} else {{
                // Odd number: middle value
                medianThreshold = featureValues[Math.floor(featureValues.length / 2)];
            }}
            
            // Bin data into low/high by feature value (using < for high to handle ties at median)
            const dataLow = allData.filter(d => d.feature < medianThreshold);
            const dataHigh = allData.filter(d => d.feature >= medianThreshold);
            
            // If all values are the same (or very skewed), split evenly to avoid empty groups
            if (dataLow.length === 0 || dataHigh.length === 0) {{
                allData.sort((a, b) => a.feature - b.feature);
                const midPoint = Math.floor(allData.length / 2);
                dataLow.length = 0;
                dataHigh.length = 0;
                dataLow.push(...allData.slice(0, midPoint));
                dataHigh.push(...allData.slice(midPoint));
            }}
            
            // Group by tool for each bin
            const tools = [...new Set(allData.map(d => d.tool))].sort();
            
            // Tool colors - kept in sync with ReportMetadata.tool_colors in create_html_report.py
            const toolColors = {{
                'MuSiCal': '#0C5DA5', 'MuSiCal-FN': '#4898CC', 'MuSiCal-BI': '#87C0E0',
                'SPA': '#E69F00', 'SPA_low': '#B36B00', 'FSPA': '#FFCC66',
                'SigID': '#007A73', 'SigTrim': '#CC3311',
                'deconstructSigs': '#8B4FCF', 'MutationalPatterns': '#5E5E5E',
                'sigfit': '#B8860B', 'yapsa': '#6A5ACD', 'sigspack': '#A0522D',
                'signature_tools': '#9932CC', 'siglasso': '#556B2F',
                'mmsig': '#C2185B', 'sigminer': '#455A64',
            }};
            
            // Create boxplot traces - separate trace for each tool, grouped by Low/High
            const traces = [];
            
            // For each tool, create a Low trace and a High trace
            tools.forEach(tool => {{
                const toolDataLow = dataLow.filter(d => d.tool === tool);
                const toolDataHigh = dataHigh.filter(d => d.tool === tool);
                
                // Low trace for this tool
                if (toolDataLow.length > 0) {{
                    traces.push({{
                        y: toolDataLow.map(d => d.metric),
                        x: Array(toolDataLow.length).fill(tool),
                        name: 'Low',
                        type: 'box',
                        boxpoints: 'all',
                        jitter: 0.3,
                        pointpos: 0,
                        marker: {{
                            size: 4,
                            opacity: 0.5,
                            color: toolColors[tool] || '#6b7280'
                        }},
                        line: {{
                            color: toolColors[tool] || '#6b7280',
                            width: 2
                        }},
                        fillcolor: `rgba${{(toolColors[tool] || '#6b7280').replace('#','').match(/.{{2}}/g).map(x=>parseInt(x,16)).join(',')}},.3)`,
                        text: toolDataLow.map(d => `Sample: ${{d.sample}}<br>${{feat}}: ${{d.feature.toFixed(3)}}<br>${{metric}}: ${{d.metric.toFixed(3)}}`),
                        hovertemplate: `<b>${{tool}} - Low</b><br>%{{text}}<extra></extra>`,
                        offsetgroup: 'Low',
                        legendgroup: 'Low',
                        showlegend: tool === tools[0]  // Only show legend for first tool
                    }});
                }}
                
                // High trace for this tool
                if (toolDataHigh.length > 0) {{
                    traces.push({{
                        y: toolDataHigh.map(d => d.metric),
                        x: Array(toolDataHigh.length).fill(tool),
                        name: 'High',
                        type: 'box',
                        boxpoints: 'all',
                        jitter: 0.3,
                        pointpos: 0,
                        marker: {{
                            size: 4,
                            opacity: 0.6,
                            color: toolColors[tool] || '#6b7280'
                        }},
                        line: {{
                            color: toolColors[tool] || '#6b7280',
                            width: 2
                        }},
                        fillcolor: `rgba${{(toolColors[tool] || '#6b7280').replace('#','').match(/.{{2}}/g).map(x=>parseInt(x,16)).join(',')}},.5)`,
                        text: toolDataHigh.map(d => `Sample: ${{d.sample}}<br>${{feat}}: ${{d.feature.toFixed(3)}}<br>${{metric}}: ${{d.metric.toFixed(3)}}`),
                        hovertemplate: `<b>${{tool}} - High</b><br>%{{text}}<extra></extra>`,
                        offsetgroup: 'High',
                        legendgroup: 'High',
                        showlegend: tool === tools[0]  // Only show legend for first tool
                    }});
                }}
            }});

            // Get theoretical minimum TAE (if available and metric is TAE)
            let theoreticalMinTAE = null;
            if (metric === 'TAE' && allData.length > 0 && allData[0].metrics && 'theoretical_min_TAE' in allData[0].metrics) {{
                theoreticalMinTAE = allData[0].metrics['theoretical_min_TAE'];
            }}
            
            // Add theoretical maximum line if applicable
            if (theoreticalMinTAE !== null && theoreticalMinTAE > 0) {{
                traces.push({{
                    x: tools,
                    y: Array(tools.length).fill(theoreticalMinTAE),
                    name: 'Theoretical Best (Unknown Sigs)',
                    type: 'scatter',
                    mode: 'lines',
                    line: {{
                        color: '#ef4444',
                        width: 2,
                        dash: 'dash'
                    }},
                    showlegend: true,
                    hovertemplate: `Theoretical Minimum TAE: ${{theoreticalMinTAE.toFixed(4)}}<br>(Due to unknown signatures)<extra></extra>`
                }});
            }}

            const layout = {{ 
                title: {{
                    text: `<b>${{metric}}</b> distribution by <b>${{feat}}</b> (Low vs High)`,
                    font: {{ size: 16, color: '#111827' }}
                }},
                xaxis: {{
                    title: {{ text: 'Tool', font: {{ size: 13, color: '#374151' }} }},
                    gridcolor: '#e5e7eb',
                    categoryorder: 'array',
                    categoryarray: tools
                }}, 
                yaxis: {{
                    title: {{ text: metric, font: {{ size: 14, color: '#374151' }} }},
                    gridcolor: '#e5e7eb',
                    zeroline: true,
                    zerolinecolor: '#9ca3af',
                    zerolinewidth: 1
                }},
                showlegend: true,
                legend: {{
                    x: 1.02,
                    y: 1,
                    bgcolor: 'rgba(255,255,255,0.9)',
                    bordercolor: '#e5e7eb',
                    borderwidth: 1,
                    title: {{ text: 'Group', font: {{ size: 12, color: '#374151' }} }}
                }},
                plot_bgcolor: 'white',
                paper_bgcolor: '#fafafa',
                margin: {{ l: 70, r: 120, t: 70, b: 90 }},
                boxmode: 'group',
                boxgap: 0.2,
                boxgroupgap: 0.1
            }};
            
            Plotly.newPlot('{widget_id}_plot', traces, layout, {{responsive:true, displayModeBar: true}});
            
            // Update stats
            document.getElementById('{widget_id}_stats').style.display = 'block';
            document.getElementById('{widget_id}_threshold').textContent = medianThreshold.toFixed(3);
            document.getElementById('{widget_id}_n_low').textContent = dataLow.length;
            document.getElementById('{widget_id}_n_high').textContent = dataHigh.length;
        }}

        document.getElementById('{widget_id}_update').addEventListener('click', updatePlot);

        // Pre-select first options and auto-plot
        if (featureSelect.options.length > 1) featureSelect.selectedIndex = 1;
        if (metricSelect.options.length > 1) metricSelect.selectedIndex = 1;
        if (mutsSelect.options.length) mutsSelect.selectedIndex = 0;
        
        // Auto-update on load
        setTimeout(updatePlot, 100);
}})();
</script>
"""
        # Prepare data for embedding (avoid CORS issues with file:// protocol)
        import json
        import pandas as pd
        
        if data_df is not None:
            # Use provided dataframe
            df = data_df
        else:
            # Fallback to reading file (but this won't work with file:// protocol)
            df = pd.read_csv(per_sample_metrics_tsv, sep='\t')
        
        # Convert dataframe to {{header: [...], rows: [...]}} format for JavaScript
        data_dict = {
            'header': df.columns.tolist(),
            'rows': df.astype(str).values.tolist()  # Convert all to strings for JSON
        }
        
        # Embed data directly in JavaScript to avoid fetch/CORS issues
        data_json = json.dumps(data_dict)
        html = html.replace('__DATA__', data_json)
        
        return html

