"""Configuration utilities for SigFitTest workflow."""

import yaml
from pathlib import Path


def get_tool_format_group(tool_name, config):
    """Get the format group information for a tool."""
    tool_formats = config.get("tool_formats", {})
    
    for format_type, format_info in tool_formats.items():
        if tool_name.lower() in [t.lower() for t in format_info.get("tools", [])]:
            return format_type, format_info
    
    return None, None


def get_tool_config(tool_name, config):
    """Get consolidated configuration for a tool."""
    # Start with defaults
    defaults = config.get("defaults", {})
    tool_config = defaults.copy()
    
    # Get format group configuration
    format_type, format_info = get_tool_format_group(tool_name, config)
    if format_info:
        tool_config.update({
            "cosmic_signatures": {
                "standard": format_info.get("standard"),
                "relevant": format_info.get("relevant")
            },
            "input_data": format_info.get("input_data_base")
        })
    
    # Apply tool-specific overrides
    tool_overrides = config.get("tool_overrides", {}).get(tool_name, {})
    if tool_overrides:
        # Deep merge for nested configs like cosmic_signatures
        for key, value in tool_overrides.items():
            if isinstance(value, dict) and key in tool_config:
                tool_config[key].update(value)
            else:
                tool_config[key] = value
    
    return tool_config


def list_tools_by_format():
    """List all tools grouped by their format requirements."""
    config_path = Path(__file__).parent.parent / "config" / "config.yaml"
    with open(config_path) as f:
        config = yaml.safe_load(f)
    
    tool_formats = config.get("tool_formats", {})
    
    print("Tools grouped by format:")
    for format_type, format_info in tool_formats.items():
        tools = format_info.get("tools", [])
        print(f"  {format_type}: {', '.join(tools)}")


if __name__ == "__main__":
    list_tools_by_format()