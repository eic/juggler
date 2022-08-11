namespace Jug {

  class Algorithm {
  private:
    std::string m_name;

    static std::function<std::ostream&()> m_info;
    static std::function<std::ostream&()> m_error;

    static void SetInfo(std::function<std::ostream&()> info) {
      m_info = info;
    }
    static void SetError(std::function<std::ostream&()> error) {
      m_error = error;
    }

  public:
    Algorithm(const std::string& name)
    : m_name ( name ) {
    }

  protected:
    static std::ostream& info() { return m_info(); };
    static std::ostream& error() { return m_error(); };

  };

} // namespace Jug
