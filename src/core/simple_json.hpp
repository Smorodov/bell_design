#pragma once

#include <cctype>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <fstream>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Small dependency-free JSON reader for the project's fixed configuration files.
// It intentionally supports the complete JSON value grammar, but exposes only
// the accessors needed by BellOptimizer. Writing is handled by the few callers
// that produce stable, human-readable result files.
namespace bell::json {

class Value {
public:
    enum class Type { Null, Bool, Number, String, Array, Object };

    Value() = default;
    explicit Value(bool v) : type_(Type::Bool), bool_(v) {}
    explicit Value(double v) : type_(Type::Number), number_(v) {}
    explicit Value(std::string v) : type_(Type::String), string_(std::move(v)) {}
    explicit Value(std::vector<Value> v) : type_(Type::Array), array_(std::move(v)) {}
    explicit Value(std::map<std::string, Value> v) : type_(Type::Object), object_(std::move(v)) {}

    Type type() const noexcept { return type_; }
    bool is_array() const noexcept { return type_ == Type::Array; }
    bool is_object() const noexcept { return type_ == Type::Object; }

    const Value* find(const std::string& key) const noexcept {
        if (!is_object()) return nullptr;
        const auto it = object_.find(key);
        return it == object_.end() ? nullptr : &it->second;
    }

    double number(double fallback) const noexcept {
        return type_ == Type::Number ? number_ : fallback;
    }

    int integer(int fallback) const noexcept {
        if (type_ != Type::Number || !std::isfinite(number_) ||
            number_ < static_cast<double>(std::numeric_limits<int>::min()) ||
            number_ > static_cast<double>(std::numeric_limits<int>::max())) return fallback;
        return static_cast<int>(number_);
    }

    bool boolean(bool fallback) const noexcept {
        return type_ == Type::Bool ? bool_ : fallback;
    }

    const std::string& string(const std::string& fallback) const noexcept {
        return type_ == Type::String ? string_ : fallback;
    }

    const std::vector<Value>& array() const noexcept { return array_; }

private:
    Type type_ = Type::Null;
    bool bool_ = false;
    double number_ = 0.0;
    std::string string_;
    std::vector<Value> array_;
    std::map<std::string, Value> object_;
};

class Parser {
public:
    explicit Parser(std::string text) : text_(std::move(text)) {}

    Value parse() {
        skip_ws();
        Value result = parse_value();
        skip_ws();
        if (pos_ != text_.size()) error("unexpected characters after JSON value");
        return result;
    }

private:
    std::string text_;
    std::size_t pos_ = 0;

    [[noreturn]] void error(const char* message) const {
        throw std::runtime_error(std::string("JSON error at position ") +
                                 std::to_string(pos_) + ": " + message);
    }

    void skip_ws() {
        while (pos_ < text_.size() && std::isspace(static_cast<unsigned char>(text_[pos_]))) ++pos_;
    }

    char peek() const {
        return pos_ < text_.size() ? text_[pos_] : '\0';
    }

    char take() {
        if (pos_ >= text_.size()) error("unexpected end of input");
        return text_[pos_++];
    }

    void expect(char c) {
        if (take() != c) error("unexpected character");
    }

    Value parse_value(std::size_t depth = 0) {
        if (depth > 256) error("maximum nesting depth exceeded");
        skip_ws();
        switch (peek()) {
            case '{': return parse_object(depth);
            case '[': return parse_array(depth);
            case '"': return Value(parse_string());
            case 't': parse_literal("true"); return Value(true);
            case 'f': parse_literal("false"); return Value(false);
            case 'n': parse_literal("null"); return Value();
            default:
                if (peek() == '-' || (peek() >= '0' && peek() <= '9')) return Value(parse_number());
                error("expected JSON value");
        }
    }

    void parse_literal(const char* literal) {
        while (*literal) {
            if (take() != *literal++) error("invalid literal");
        }
    }

    Value parse_object(std::size_t depth) {
        expect('{');
        skip_ws();
        std::map<std::string, Value> object;
        if (peek() == '}') { ++pos_; return Value(std::move(object)); }
        for (;;) {
            skip_ws();
            if (peek() != '"') error("object key must be a string");
            std::string key = parse_string();
            skip_ws();
            expect(':');
            object.emplace(std::move(key), parse_value(depth + 1));
            skip_ws();
            const char c = take();
            if (c == '}') break;
            if (c != ',') error("expected ',' or '}'");
        }
        return Value(std::move(object));
    }

    Value parse_array(std::size_t depth) {
        expect('[');
        skip_ws();
        std::vector<Value> array;
        if (peek() == ']') { ++pos_; return Value(std::move(array)); }
        for (;;) {
            array.push_back(parse_value(depth + 1));
            skip_ws();
            const char c = take();
            if (c == ']') break;
            if (c != ',') error("expected ',' or ']'");
        }
        return Value(std::move(array));
    }

    std::string parse_string() {
        expect('"');
        std::string result;
        while (pos_ < text_.size()) {
            const char c = take();
            if (c == '"') return result;
            if (static_cast<unsigned char>(c) < 0x20) error("control character in string");
            if (c != '\\') { result += c; continue; }

            const char escape = take();
            switch (escape) {
                case '"': result += '"'; break;
                case '\\': result += '\\'; break;
                case '/': result += '/'; break;
                case 'b': result += '\b'; break;
                case 'f': result += '\f'; break;
                case 'n': result += '\n'; break;
                case 'r': result += '\r'; break;
                case 't': result += '\t'; break;
                case 'u': parse_unicode_escape(result); break;
                default: error("invalid string escape");
            }
        }
        error("unterminated string");
    }

    void parse_unicode_escape(std::string& out) {
        unsigned code = 0;
        for (int i = 0; i < 4; ++i) {
            const char c = take();
            code <<= 4;
            if (c >= '0' && c <= '9') code += static_cast<unsigned>(c - '0');
            else if (c >= 'a' && c <= 'f') code += static_cast<unsigned>(c - 'a' + 10);
            else if (c >= 'A' && c <= 'F') code += static_cast<unsigned>(c - 'A' + 10);
            else error("invalid unicode escape");
        }
        if (code <= 0x7f) out += static_cast<char>(code);
        else if (code <= 0x7ff) {
            out += static_cast<char>(0xc0 | (code >> 6));
            out += static_cast<char>(0x80 | (code & 0x3f));
        } else {
            out += static_cast<char>(0xe0 | (code >> 12));
            out += static_cast<char>(0x80 | ((code >> 6) & 0x3f));
            out += static_cast<char>(0x80 | (code & 0x3f));
        }
    }

    double parse_number() {
        const std::size_t start = pos_;
        if (peek() == '-') ++pos_;
        if (peek() == '0') ++pos_;
        else {
            if (peek() < '1' || peek() > '9') error("invalid number");
            while (peek() >= '0' && peek() <= '9') ++pos_;
        }
        if (peek() == '.') {
            ++pos_;
            if (peek() < '0' || peek() > '9') error("invalid number fraction");
            while (peek() >= '0' && peek() <= '9') ++pos_;
        }
        if (peek() == 'e' || peek() == 'E') {
            ++pos_;
            if (peek() == '+' || peek() == '-') ++pos_;
            if (peek() < '0' || peek() > '9') error("invalid number exponent");
            while (peek() >= '0' && peek() <= '9') ++pos_;
        }
        char* end = nullptr;
        const double value = std::strtod(text_.c_str() + start, &end);
        if (end != text_.c_str() + pos_) error("invalid number");
        if (!std::isfinite(value)) error("number out of representable range");
        return value;
    }
};

inline Value parse_file(const std::string& filename) {
    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open()) throw std::runtime_error("cannot open file: " + filename);
    return Parser(std::string((std::istreambuf_iterator<char>(file)),
                              std::istreambuf_iterator<char>())).parse();
}

inline double number(const Value& object, const char* key, double fallback) {
    const Value* value = object.find(key);
    return value ? value->number(fallback) : fallback;
}

inline int integer(const Value& object, const char* key, int fallback) {
    const Value* value = object.find(key);
    return value ? value->integer(fallback) : fallback;
}

inline bool boolean(const Value& object, const char* key, bool fallback) {
    const Value* value = object.find(key);
    return value ? value->boolean(fallback) : fallback;
}

inline std::string string(const Value& object, const char* key, std::string fallback) {
    const Value* value = object.find(key);
    return value ? value->string(fallback) : fallback;
}

inline void number_array(const Value& object, const char* key, double* out, std::size_t count) {
    const Value* value = object.find(key);
    if (!value || !value->is_array()) return;
    const auto& array = value->array();
    const std::size_t n = array.size() < count ? array.size() : count;
    for (std::size_t i = 0; i < n; ++i) out[i] = array[i].number(out[i]);
}

} // namespace bell::json
